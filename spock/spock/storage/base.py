from abc import ABC, abstractmethod
from typing import Literal, Iterable, Iterator

from rdkit import Chem

from spock.utils.standardizers.base import ChemStandardizer
from spock.docking.base import Protein


class StoredMol(ABC):
    """
    A simple interface for a molecule that can be stored in a chemstore.
    """

    @property
    @abstractmethod
    def smiles(self):
        """
        Get the SMILES of the molecule.

        :return: SMILES of the molecule
        """
        pass

    @property
    @abstractmethod
    def id(self):
        """
        Get the identifier of the molecule.

        :return: identifier of the molecule
        """
        pass

    @property
    @abstractmethod
    def metadata(self):
        """
        Get the metadata of the molecule.

        :return: metadata of the molecule
        """
        pass

    @metadata.setter
    @abstractmethod
    def metadata(self, value):
        """
        Set the metadata of the molecule.

        :param value: metadata of the molecule
        """
        pass

    @property
    @abstractmethod
    def representations(self) -> list[Chem.Mol]:
        pass

    def as_rd_mol(self) -> Chem.Mol:
        """
        Get the RDKit molecule object of the standardized representation of this instance.

        :return: `rdkit.Chem.Mol` instance
        """

        return Chem.MolFromSmiles(self.smiles)


class ChemIdentifier(ABC):

    @abstractmethod
    def __call__(self, smiles: str) -> str:
        """
        Get the identifier of the molecule represented by the given SMILES.

        :param smiles: input SMILES
        :return: calculated identifier
        """
        pass


class InchiIdentifier(ChemIdentifier):

    def __call__(self, smiles: str) -> str:
        """
        Get the InChIKey of the molecule represented by the given SMILES.

        :param smiles: input SMILES
        :return: calculated InChIKey
        """
        return Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))


class ChemStore(ABC):

    @property
    @abstractmethod
    def standardizer(self) -> ChemStandardizer:
        """
        Get the standardizer used by the store.

        :return: `ChemStandardizer` instance
        """
        pass

    @property
    @abstractmethod
    def identifier(self) -> ChemIdentifier:
        """
        Get the identifier used by the store.

        :return: `ChemIdentifier` instance
        """
        pass

    @property
    def n_mols(self) -> int:
        """Number of molecules in storage."""
        return self.get_mol_count()

    @abstractmethod
    def get_mol(self, mol_id: str) -> StoredMol:
        """
        Get a molecule from the store using its ID.

        :param mol_id: identifier of the molecule to search
        :return: instance of `StoredMol`
        """
        pass

    def add_mol_from_smiles(
        self,
        smiles,
        metadata: dict = None,
        library=None,
        raise_on_existing=False,
        update_existing: Literal["append", "update", None] = None,
        **kwargs,
    ):
        """
        Add a molecule to the store using its raw SMILES. The SMILES will be standardized and an identifier will be
        calculated.

        :param smiles: SMILES of the molecule to add
        :param metadata: additional metadata to store with the molecule
        :param library: library to add the molecule to
        :param raise_on_existing: raise an exception if the molecule already exists in the store
        :param update_existing: update the metadata of the existing molecule if it already exists in the store (if `append` is specified, the metadata will be appended to the existing metadata, if `update` is specified, the existing metadata will be overwritten)
        :param kwargs: additional arguments to pass to `self.add_mol` after standardization and identification

        :return: `StoredMol` instance of the added molecule

        :raises: `ValueError` if the molecule cannot be added or if `raise_on_existing` is set to `True` and the molecule already exists in the store
        """
        metadata = metadata or dict()
        metadata["standardizer_settings"] = self.standardizer.settings
        metadata["original_smiles"] = smiles
        smiles, _ = self.standardizer.convert_smiles(smiles)
        identifier = self.identifier(smiles)
        return self.add_mol(
            smiles,
            identifier,
            metadata,
            library=library,
            raise_on_existing=raise_on_existing,
            update_existing=update_existing,
            **kwargs,
        )

    @abstractmethod
    def add_mol(self, smiles, mol_id, *args, **kwargs):
        """
        Add a molecule to the store. This method should not perform any standardization or identifier calculation. The `add_mol_from_smiles` method should be used instead if automatic standardization and identification should be performed before storage.

        :param smiles: molecule to add as SMILES
        :param mol_id: identifier of the molecule to add as determined by `self.identifier`
        :param metadata: additional metadata to store with the molecule
        :return: `StoredMol` instance of the added molecule

        :raises: `ValueError` if the molecule cannot be added
        """
        pass

    @abstractmethod
    def remove_mol(self, mol_id):
        """
        Remove a molecule from the store.

        :param mol_id: identifier of the molecule to remove
        :return:
        """
        pass

    @abstractmethod
    def get_mol_ids(self):
        """
        Get all molecule IDs in the store.

        :return: list of molecule IDs
        """
        pass

    @abstractmethod
    def get_mol_count(self):
        """
        Get the number of molecules in the store.

        :return: number of molecules
        """
        pass

    @abstractmethod
    def iter_mols(self) -> Iterator:
        """
        Iterate over all molecules in the store.

        :return: iterator over `StoredMol` instances
        """
        pass

    def iter_chunks(self, size=1000) -> Iterator[list[StoredMol]]:
        """
        Iterate over chunks of molecules across the store.

        :return: an iterable of lists of stored molecules
        """
        ret = []
        for mol in self.iter_mols():
            ret.append(mol)
            if len(ret) % size == 0:
                yield ret
                ret = []
        yield ret

    @abstractmethod
    def add_sdf(self, sdf_path):
        pass

    @abstractmethod
    def save(self, *args, **kwargs):
        pass

    @abstractmethod
    def clear_sdfs(self):
        pass

    @abstractmethod
    def clear_poses(self, target: Protein = None):
        pass

    @abstractmethod
    def add_target(self, target: Protein):
        """
        Add a target for the molecule.

        :param target: target to add
        """
        pass

    @property
    @abstractmethod
    def targets(self):
        pass

    @abstractmethod
    def add_poses(self, mol_id: str, repr_index: int, poses: Chem.Mol, target: Protein):
        """
        Add docking poses to the given molecule.

        :param mol_id: id of the molecule to add the poses to
        :param repr_index: index of the representation the poses were generated with (see `StoredMol.representations`)
        :param poses: poses to add (this is an RDKit molecule with the docking poses stored in the conformers)
        :param target: the target to which the poses correspond to
        """
        pass

    @abstractmethod
    def get_poses(self, mol_id: str, target: Protein):
        pass

    @abstractmethod
    def get_target(self, name: str):
        pass

    def __len__(self):
        return self.get_mol_count()

    def __contains__(self, item):
        return item in self.get_mol_ids()

    def __iter__(self):
        return self.iter_mols()

    def __getitem__(self, item):
        return self.get_mol(item)

    def __delitem__(self, key):
        return self.remove_mol(key)

    def __setitem__(self, key, value):
        return self.add_mol(value, key)

    def __str__(self):
        return f"{self.__class__.__name__} ({self.get_mol_count()})"

    def __repr__(self):
        return str(self)

    def __bool__(self):
        return len(self) > 0


class SMARTSSearchable(ABC):

    @abstractmethod
    def smarts_search(self, smarts: str):
        pass
