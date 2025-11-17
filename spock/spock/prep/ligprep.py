import math
import os
import warnings

from tqdm import tqdm

from spock.config import LigprepParams
from spock.prep.base import LigandPreparation
from spock.storage.base import ChemStore
from spock.utils.exec import run_executable

# TODO:
# 1. Fix the workers issue - overbooks licenses


class Ligprep(LigandPreparation):
    def __init__(self, schrodinger_path: str, store: ChemStore=None, config: LigprepParams=None,
                 config_path: str=None):

        self.executable = os.path.join(schrodinger_path, "ligprep")
        self.params = config if config is not None else LigprepParams(path=config_path, schrodinger_path=schrodinger_path, executable=f"{schrodinger_path}/ligprep")

    def help(self):
        """
        Print the help message for ligprep
        """
        command = f"{self.executable} -long_help"
        stdout, stderr = run_executable(command)
        print(stdout.decode("utf-8"), stderr.decode("utf-8"))

    def run(self, filename, out_path=None, output_format="sd", wait=False):
        """
        Run ligprep on a single molecule
        """
        filename = os.path.abspath(filename)
        input_extension = filename.split(".")[-1]
        in_path = os.path.dirname(filename)
        in_file = os.path.basename(filename).split(".")[0]
        output_filename = filename.split("/")[-1]
        output_filename = f"{in_file}.{output_format}"
        formats = ["smi", "sdf", "mol2", "pdb", "mol", "mae", "maegz", "mae.gz", "csv"]

        if input_extension not in formats:
            raise ValueError(f"Unsupported input format: {input_extension}")

        ligprep_kwargs = dict(
            input_string = f"-i{input_extension} {filename}",
            output_format = f"-o{output_format} {output_filename}",
            num_stereoisomers = f"-s {self.params.num_stereoisomers}",
            max_atoms = f"-m {self.params.max_atoms}",
            force_field = f"-bff {self.params.force_field}",
            epik = f"-{self.params.epik}" if self.params.epik else "",
                        optimize = "-cgx_noopt" if not self.params.optimize else "",
                        ph = f"-ph {self.params.ph}",
                                i= f"-i {self.params.i}",
            HOST = f"-HOST localhost:{self.params.cpus}",
            NJOBS = f"-NJOBS {self.params.workers}",
            wait = "-WAIT" if wait else "",
        )

        command = f'{self.executable} {" ".join(ligprep_kwargs.values())}'
        stdout, stderr = run_executable(command, cwd=in_path)
        if stderr != b"":
            print(stdout.decode("utf-8"), stderr.decode("utf-8"))

        # Move the output file to the output folder
        if out_path is None:
            # out_path = os.path.join(os.path.dirname(in_path), 'sdf')
            out_path = os.path.join(in_path, in_file+".sd")
        ########################
        # out_path = os.path.join(out_path, output_filename)
        # ligprep_output_file = os.path.join(in_path, output_filename)
        # if os.path.exists(ligprep_output_file):
        #     shutil.move(ligprep_output_file, out_path)
        # else:
        #     out_path = None
        ########################
        return out_path

    def process_store(self, store, cpus=None, chunks=20, workers=1, **ligprep_kwargs):
        """
        Process all molecules in a store
        """
        new_df = store._df.copy()
        new_df.drop(columns=["metadata", "original_smiles", "library", "sdf_count"], inplace=True)
        ## Make 'SMILES' the first column
        cols = new_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        new_df = new_df[cols]
        new_df["NAME"] = new_df["id"]

        self.params.workers = workers
        self.params.cpus = cpus if cpus is not None else self.params.cpus

        print("Running ligprep with the following parameters on the store:")
        print(self.params.as_dict())

        # Split the csv file into chunks
        mols_per_chunk = math.ceil(len(new_df)/chunks)
        for i in tqdm(range(chunks), desc="Processing files...", total=chunks):
            chunk = new_df.iloc[i*mols_per_chunk:(i+1)*mols_per_chunk]
            prefix = "store_ligprep"
            path = os.path.join(store.path, f"{prefix}.csv")
            chunk[["SMILES","NAME", "id"]].to_csv(path, index=False, header=True)
            output_sd_path = self.run(path, wait=True)
            err_path = path.replace(".csv", f"_failed_{i}.csv")
            if os.path.exists(f"{store.path}/{prefix}-dropped-indices.txt"):
                warnings.warn(f"Failed to ligprep chunk: {i}")
                os.rename(path, err_path)
                os.rename(f"{store.path}/{prefix}.log", f"{store.path}/{prefix}_failed_{i}.log")
                # with open(f"{store.path}/{prefix}-dropped-indices.txt", "r") as f:
                #     indices_failed = f.read().strip().replace("smiconvert_in:", "").strip().split(",")
                #     for idx in indices_failed:
                #         idx = int(idx.strip())
                #         warnings.warn(f"Failed to ligprep molecule: {chunk.iloc[idx,:].NAME} ({chunk.iloc[idx,:].SMILES}) in chunk={i}, see {store.path}/{prefix}_{i}.log", stacklevel=2)
                os.rename(f"{store.path}/{prefix}-dropped.csv", f"{store.path}/{prefix}-dropped_failed_{i}.csv")
                os.rename(f"{store.path}/{prefix}-dropped-indices.txt", f"{store.path}/{prefix}-dropped-indices_failed_{i}.txt")
            if not os.path.exists(output_sd_path):
                if not os.path.exists(err_path):
                    os.rename(path, err_path)
                warnings.warn(f"Output file from ligprep not generated  for chunk={i}. Bad molecules saved in: {err_path}")
                continue
            store.add_sdf(output_sd_path)
            os.remove(output_sd_path)
            if os.path.exists(f"{store.path}/{prefix}.log"):
                os.remove(f"{store.path}/{prefix}.log")
            if os.path.exists(f"{store.path}/{prefix}.csv"):
                os.remove(f"{store.path}/{prefix}.csv")


    def split_sdf(self, sdf_path, output_folder):
        """
        Split an sdf file into individual sdf files
        """

    def _process_stored_molecule(self, mol, output_format="sd"):
        """
        Process a single molecule stored in the database
        """
        path = mol.metadata["sdf_count"]
        if os.path.exists(path):
            return path

        # Generate the sdf file
        csv_path = os.path.join(mol.store_path, "smiles", mol.id + ".csv")
        if not os.path.exists(csv_path):
            mol.to_file(os.path.join(mol.store_path, "smiles"))

        # Run ligprep
        sdf_path = self.run(csv_path,
                            os.path.join(mol.store_path, "sdf"),
                            self.params.force_field,
                            self.params.max_atoms,
                            self.params.epik,
                            self.params.delete,
                            self.params.num_stereoisomers)
        return sdf_path







### OLD CODE ###
    # def process_store(self, store, chunk_size=32, sub_chunk=16,
    #                   workers=8, **ligprep_kwargs):
    #     """
    #     Process all molecules in a store

    #     Arguments:
    #         store (Store)           : store containing the molecules
    #         chunk_size (int)        : number of molecules to process in a single chunk
    #         sub_chunk (int)         : number of chunks to process in a single sub_chunk
    #         workers (int)           : number of parallel workers
    #         **ligprep_kwargs        : ligprep arguments
    #     Returns:
    #         None                   To be implemented
    #     """

    #     self.configure(**ligprep_kwargs)

    #     ids = range(chunk_size*sub_chunk)
    #     filenames = [os.path.join(store.path, 'smiles', f'{id}.csv') for id in ids]

    #     # Write the SMILES to files (partitioned for the parallel workers)
    #     for i, phile in enumerate(filenames):
    #         num_ligands = math.ceil(len(store._df) / (chunk_size*sub_chunk))
    #         idx = slice(i*num_ligands, (i+1)*num_ligands)
    #         smiles = store._df.SMILES.values[idx]
    #         ids = store._df.id.values[idx]
    #         with open(phile, 'w') as f:
    #             f.write('SMILES,id\n')
    #             for smile, id in zip(smiles, ids):
    #                 f.write(f'{smile},{id}\n')

    #     # Partition into chunks for the progress bar
    #     subdivisions = [filenames[i*chunk_size:(i+1)*chunk_size] for i in range(sub_chunk)]
    #     pool = Pool
    #     try:
    #         with Pool(workers) as pool:
    #             for _, chunk in tqdm(enumerate(subdivisions), total=sub_chunk, unit='chunk'):
    #                 _ = pool.map(self.run, chunk)
    #     except KeyboardInterrupt:
    #         print('Interrupted by user.')

    #     pool.terminate()
    #     pool.join()

    #     ## Add to the store


    #     return None
