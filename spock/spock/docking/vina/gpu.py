import os, time, datetime
import docker
import random
from spock.docking.vina.base import BaseVinaRunner
from spock.docking.vina.utils import process_stdout
from spock.utils.logs.logs import write_to_log


class VinaGPU(BaseVinaRunner):
    """
    Class methods for running Vina-GPU docker container
    Also contains methods for preparing the ligand and target:
        - Ligand preparation via rdkit and meeko
        - Target preparation via ADFR Suite and pdb_tools
    """
    def __init__(self, docker_image_name='andriusbern/vina-gpu', devices=['0'], visualize=False):

        super(VinaGPU, self).__init__(device='gpu')

        self.visualize = visualize
        self.device_id = devices

        ## Configuration for running the Vina-GPU docker container 
        # (requires nvidia-docker runtime)
        self.container = None
        dev_req = docker.types.DeviceRequest  # type: ignore
        self.docker_kwargs = dict(
            image=docker_image_name,
            runtime='nvidia',    # Use nvidia-docker runtime
            volumes = [f'{self.out_path}:{self.docking_dir}'],
            device_requests=[dev_req(device_ids=devices, capabilities=[['gpu']])])
        

    def dock(self, target_pdb_path, smiles=[], ligand_pdbqt_paths=[], output_subfolder='', 
             box_center=(0,0,0), box_size=(20,20,20), search_depth=3,
             threads=256, threads_per_call=256, clean=True, verbose=True, 
             visualize_in_pymol=False, write_log=True, **kwargs):
        """
        Use Vina-GPU docker image to dock ligands (list of SMILES or .pdbqt files) to the target. 
        Produces a .pdbqt file for each ligand (with multiple docked orientations). 

        Arguments:
            target_pdb_path (str)                   : path to target pdb file
            smiles: (list(str))                     : list of smiles strings    
            ligand_pdbqt_paths (list(str))          : list of paths to ligand pdbqt files
            output_subfolder (str), opt             : subfolder to save output files
            active_site_coords (tuple(float)), opt  : coordinates of the active site of the target (x,y,z)=(0,0,0)
            bbox_size (tuple(float)), opt           : size of the bounding box around the active site (x,y,z)=(20,20,20)
            threads (int), opt                      : number of threads to use for docking
            thread_per_call (int), opt              : number of threads to use for each call to Vina
            clean (bool), opt                       : remove ligand .pdbqt files after docking
            verbose (bool), opt                     : print docking progress, scores, etc.
            visualize_in_pymol (bool), opt          : visualize the docking results in pymol
            write_log (bool), opt                   : write log file with docking results
        Returns:
            all_scores (list(list((float)))         : list of docking scores for each ligand
        """

        assert (len(ligand_pdbqt_paths) > 0) or (len(smiles) > 0), \
        "Either a list of ligand .pdbqt paths or a list of smiles strings must be provided"

        results_path = os.path.join(self.out_path, output_subfolder)
        os.makedirs(results_path, exist_ok=True)

        # Prepare target .pdbqt file
        target_pdbqt_path = self.prepare_target(target_pdb_path, output_path=results_path)

        # Prepare ligand .pdbqt files
        print('Processing ligands...') if verbose else None
        for i, mol in enumerate(smiles): 
            uid = random.randint(0, 1000000)   
            ligand_pdbqt_path = os.path.join(results_path, f'ligand_{uid}.pdbqt')
            out_path = self.prepare_ligand(mol, out_path=ligand_pdbqt_path)
            if out_path is not None:
                ligand_pdbqt_paths.append(ligand_pdbqt_path)
        basenames = [os.path.basename(p) for p in ligand_pdbqt_paths] # Ligand basenames (format 'ligand_0.pdbqt')
        basenames_docked = [lig.replace('.pdbqt', '_docked.pdbqt') for lig in basenames] # Docked ligand basenames (format 'ligand_0_docked.pdbqt')
        ligand_paths_docked = [os.path.join(results_path, p) for p in basenames_docked]
        
        ### Start Vina-GPU docker container
        self.container = self.start_docker_container()
        try:
            timing, dates = [], []
            all_scores = [[0] for i in range(len(smiles))]
            target = os.path.basename(target_pdb_path).strip('.pdbqt')
            for i, ligand_file in enumerate(basenames):
                t0 = time.time()

                docking_args = dict(
                    receptor = f'docking/{output_subfolder}/{os.path.basename(target_pdbqt_path)}',
                    ligand   = f'docking/{output_subfolder}/{ligand_file}',
                    out      = f'docking/{output_subfolder}/{basenames_docked[i]}',
                    center_x = box_center[0],
                    center_y = box_center[1],
                    center_z = box_center[2],
                    size_x   = box_size[0],
                    size_y   = box_size[1],
                    size_z   = box_size[2],
                    thread   = threads,
                    search_depth = search_depth,
                    thread_per_call = threads_per_call)

                cmd = './Vina-GPU ' + ' '.join([f'--{k} {v}' for k, v in docking_args.items()])

                try:
                    _, (stdout, stderr) = self.container.exec_run(
                        cmd=cmd,
                        workdir=self.vina_dir,
                        demux=True)

                    scores = process_stdout(stdout)

                    if len(scores) > 0 and scores != [None]:
                        all_scores[i] = scores

                    timing += [round(time.time() - t0, 2)]
                    dates += [datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")]
                    if verbose:
                        print(f'- {self.device}:{self.device_id} | [{dates[-1]} | t={timing[-1]}s] Docked ligand {i+1}/{len(basenames)} | Affinity values: {all_scores[i]}...')
                    
                    if write_log:
                        log_path = os.path.join(results_path, 'log.tsv')
                        write_to_log(log_path, smiles[i], target, all_scores[i], ligand_paths_docked[i])

                    if clean: # Remove intermediate files (undocked ligand .pdbqt files)
                        os.remove(ligand_pdbqt_paths[i])
                        os.remove(ligand_paths_docked[i])
                except Exception as d:
                    print(d)
                    
        except Exception as e:
            print(f'Error has occurred while docking ligand {i}: {e, stderr}')
            raise e
        except KeyboardInterrupt:
            print('Docking interrupted by user')
        finally:
            self.remove_docker_container()

        if visualize_in_pymol or self.visualize: 
            self.visualize_results(target_pdb_path, ligand_paths_docked, scores=all_scores)
    
        return all_scores
