"""
Contains the docking base class
"""

import logging
import os
import time
from abc import ABC, abstractmethod
from rdkit import Chem

from spock.parallel.iterators import parallel_iter, batched_iter


class Protein(ABC):
    """
    A protein object
    """

    def __init__(self, name: str, *args, **kwargs):
        self.name = name

    @abstractmethod
    def get_pdb(self):
        pass

    def get(self):
        return Chem.MolFromPDBBlock(self.get_pdb())


class ProteinFromPDBString(Protein):
    """
    A protein object from a PDB string
    """

    def __init__(self, name: str, pdb_string: str, *args, **kwargs):
        """
        Creates a protein object.

        :param pdb_string: The PDB string.
        """
        super().__init__(name, *args, **kwargs)
        self.pdb_string = pdb_string

    def get_pdb(self):
        return self.pdb_string


class ProteinFromPDBFile(Protein):
    """
    A protein object from a PDB file
    """

    def __init__(self, name: str, pdb_file: str, *args, **kwargs):
        """
        Creates a protein object.

        :param pdb_file: Path to the file.
        """
        super().__init__(name, *args, **kwargs)
        self.pdb_file = pdb_file

    def get_pdb(self, raw=True):
        with open(self.pdb_file, "r") as f:
            pdb = f.read()
            return pdb


class Docking(ABC):

    def __init__(self, protein: Protein, *args, **kwargs):
        self.protein = protein

    @abstractmethod
    def dock_mol(self, mol: Chem.Mol) -> tuple[Chem.Mol, tuple[dict]]:
        pass

    def dock_stored(self, mol: "StoredMol"):
        repr_poses = []
        for mol_repr in mol.representations:
            data = {"poses": None, "metadata": None}
            poses, metadata = self.dock_mol(mol_repr)
            data["poses"] = poses
            data["metadata"] = metadata
            repr_poses.append(data)
        return repr_poses

    def dock_stored_list(self, mols: list["StoredMol"]):
        ret = dict()
        for mol in mols:
            ret[mol.id] = self.dock_stored(mol)
        return ret

    def dock_smiles_list(self, smiles: list[str]):
        return self.dock_molecule_list([Chem.MolFromSmiles(x) for x in smiles])

    def dock_molecule_list(self, mols: list[Chem.Mol]):
        poses = []
        for mol in mols:
            poses.append(self.dock_mol(mol))
        return poses

    def dock_storage(self, storage: "ChemStore", overwrite=False, save=False):
        if overwrite:
            storage.clear_poses(self.protein)
        for mol in storage:
            poses = self.dock_stored(mol)
            for idx_repr, poses_repr in enumerate(poses):
                storage.add_poses(mol.id, idx_repr, poses_repr, self.protein)
            if save:
                storage.save()


class ParallelDocking(Docking, ABC):

    def __init__(self, protein, *args, n_cpus=None, timeout=None, **kwargs):
        super().__init__(protein, *args, **kwargs)
        self.n_cpus = n_cpus if n_cpus is not None else os.cpu_count()
        self.stats = None
        self.timeout = timeout

    def _dock_mols(self, mols: list[Chem.Mol]):
        ret = []
        for mol in mols:
            ret.append(self.dock_mol(mol))
        return ret

    def dock_molecule_list(self, mols: list[Chem.Mol], chunk_size=100):
        poses = []
        for result in parallel_iter(
            batched_iter(mols, chunk_size),
            self._dock_mols,
            self.n_cpus,
            timeout=self.timeout,
        ):
            if isinstance(result, Exception):
                logging.error(result, exc_info=True)
                continue
            poses.extend(result)
        return poses

    def dock_storage(
        self, storage: "ChemStore", overwrite=False, save=False, chunk_size=1
    ):
        if overwrite:
            storage.clear_poses(self.protein)
        before_time = time.perf_counter()
        self.stats = {
            "chunk_times": [],
        }
        for result in parallel_iter(
            storage.iter_chunks(chunk_size),
            self.dock_stored_list,
            self.n_cpus,
            timeout=self.timeout,
        ):
            if isinstance(result, Exception):
                logging.error(result, exc_info=True)
                continue
            for mol_id, poses_info in result.items():
                for idx_repr, info_repr in enumerate(poses_info):
                    poses = info_repr["poses"]
                    if poses is None:
                        logging.error(
                            f"Failed to generate poses "
                            f"for molecule {mol_id} with SMILES: "
                            f"{storage.get_mol(mol_id).smiles}."
                        )
                        continue
                    for conf_id in range(poses.GetNumConformers()):
                        conf = poses.GetConformer(conf_id)
                        for key in info_repr["metadata"][conf_id]:
                            conf.SetProp(key, str(info_repr["metadata"][conf_id][key]))
                    storage.add_poses(mol_id, idx_repr, poses, self.protein)
            if save:
                storage.save()
            chunk_time = time.perf_counter() - before_time
            self.stats["chunk_times"].append(chunk_time)
            before_time = time.perf_counter()
        print(f"Docking of '{storage}' finished.")
        chunk_times = self.stats["chunk_times"]
        print(f"Average time per chunk: {sum(chunk_times) / len(chunk_times)} s")
