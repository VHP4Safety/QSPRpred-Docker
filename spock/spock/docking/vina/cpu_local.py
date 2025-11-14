import logging

from vina import Vina
from rdkit import Chem

from .utils import prep_mol, parse_poses
from ..base import ParallelDocking, ProteinFromPDBFile


class VinaProtein(ProteinFromPDBFile):

    def __init__(self, name: str, pdb_file: str, pdbqt_file: str, *args, **kwargs):
        """
        Creates a protein object.

        :param pdb_file: Path to the file.
        """
        super().__init__(name, pdb_file, *args, **kwargs)
        self.pdbqt_file = pdbqt_file


class VinaDockingCPULocal(ParallelDocking):

    def __init__(
        self,
        protein: VinaProtein,
        box_spec: dict,
        n_cpus: int = None,
        embed_mols=True,
        add_hs=True,
        mol_prep=prep_mol,
        pose_parser=parse_poses,
        timeout: int = None,
        seed: int = None,
        *args,
        **kwargs
    ):
        super().__init__(protein, *args, n_cpus=n_cpus, timeout=timeout, **kwargs)
        self.protein = protein
        self.box_spec = box_spec
        self.args = args
        self.kwargs = kwargs
        self.embed_mols = embed_mols
        self.add_hs = add_hs
        self.mol_prep = mol_prep
        self.pose_parser = pose_parser
        self.seed = seed

    def get_vina(self):
        vina = Vina(sf_name="vina", cpu=1, seed=self.seed)
        vina.set_receptor(self.protein.pdbqt_file)
        return vina

    def dock_mol(self, mol: Chem.Mol) -> tuple[Chem.Mol, tuple[dict]]:
        try:
            vina = self.get_vina()
            mol_pdbqt = self.mol_prep(
                mol, embed=self.embed_mols, add_hs=self.add_hs, seed=self.seed
            )
            vina.set_ligand_from_string(mol_pdbqt)
            vina.compute_vina_maps(
                center=self.box_spec["center"], box_size=self.box_spec["box_size"]
            )
            vina.dock(*self.args, **self.kwargs)
            poses = self.pose_parser(vina.poses())
            n_poses = poses.GetNumConformers()
            energies = vina.energies(n_poses)
            metadata = []
            for conf_id in range(poses.GetNumConformers()):
                metadata.append(
                    {
                        "vina_energy_total": str(energies[conf_id, 0]),
                        "vina_energy_inter": str(energies[conf_id, 1]),
                        "vina_energy_intra": str(energies[conf_id, 2]),
                        "vina_energy_torsions": str(energies[conf_id, 3]),
                        "vina_energy_intra_best_pose": str(energies[conf_id, 4]),
                    }
                )
            return poses, metadata
        except Exception as e:
            logging.exception(e, exc_info=True, stack_info=True)
            return None, None
