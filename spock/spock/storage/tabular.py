import io
import json
import logging
import warnings

import meeko
import pandas as pd
import os
from rdkit import Chem
from tqdm import tqdm

from spock.parallel.iterators import parallel_iter
from spock.storage.base import (
    StoredMol,
    ChemStandardizer,
    ChemStore,
    InchiIdentifier,
    SMARTSSearchable,
)
from spock.utils.molecules import smarts_match
from spock.utils.standardizers.papyrus import PapyrusStandardizer
from spock.parallel import parallel_apply
from spock.utils.formats.mol import view_difference
from spock.utils.dataframes import load_dataframe
from spock.utils.formats import sdf_iterator, parse_sdf, sdf_to_lines
from spock.utils import encode, decode
from typing import List, Any, Iterator, Literal
from spock.docking.base import Protein, ProteinFromPDBString

### TODO:
# - Add standardizer metadata back into the dataframe (as a ref to a config file)
# - Metadata update after parallel standardization
# - Add a method to update the dataframe with new columns


def read_pdb(pdb_file):
    # pdb_file = open(pdb_file, "r").read()
    # prot = Chem.MolFromMolBlock(pdb_file, sanitize=False, removeHs=False)
    prot = Chem.MolFromPDBFile(
        pdb_file, removeHs=False, sanitize=False, proximityBonding=True
    )
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SYMMRINGS)
    # Chem.SanitizeMol(prot, Chem.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(prot, Chem.SANITIZE_CLEANUP)
    # prot = Chem.AddHs(prot, addCoords=True)
    return prot


def parse_poses(poses):
    pmol = meeko.PDBQTMolecule(poses)
    return meeko.RDKitMolCreate.from_pdbqt_mol(pmol)[0]


class TabularMol(StoredMol):

    def __init__(self, parent, mol_id: str, smiles=None, index=None) -> None:
        """
        Create a new molecule instance.

        :param parent: parent ChemStore instance
        :param mol_id: identifier of the molecule
        :param smiles: SMILES of the molecule
        :param index: index of the molecule in the dataframe
        """
        self.parent = parent
        self.index = (
            parent._df[parent._df.id == mol_id].index[0] if index is None else index
        )
        self._id = mol_id
        self._smiles = (
            parent._df[parent._df.id == mol_id]["SMILES"].values[0]
            if smiles is None
            else smiles
        )

    def __repr__(self) -> str:
        return f"TabularMol({self.id}, {self.smiles})"

    @property
    def id(self) -> str:
        """
        Get the identifier of the molecule.
        """
        return self._id

    @property
    def smiles(self) -> str:
        """
        Get the SMILES of the molecule.
        """
        return self._smiles

    @property
    def metadata(self) -> pd.Series:
        """
        Get the row of the dataframe corresponding to this molecule.
        """
        return self.parent._df.iloc[self.index]

    @metadata.setter
    def metadata(self, data: pd.Series) -> None:
        """
        Set the row of the dataframe corresponding to this molecule.
        """
        self.parent._df.iloc[self.index] = data

    @property
    def representations(self):
        sdfs = self.sdf()
        ret = []
        for sdf in sdfs:
            mol = Chem.MolFromMolBlock(
                sdf, strictParsing=False, sanitize=False, removeHs=False
            )
            properties = parse_sdf(next(sdf_to_lines(sdf.split("\n"))))
            for prop in properties:
                mol.SetProp(prop, properties[prop])
            ret.append(mol)
        return ret

    def add_data(self, property: str, data: Any) -> None:
        """
        Add data to the dataframe.

        :param property: name of the column
        :param data: data to be added
        """
        self.parent._df.loc[self.index, property] = data
        # self.parent._df.iloc[self.index][property] = data

    def get_data(self, property: str) -> Any:
        """
        Get data from the dataframe.

        :param property: name of the column
        :return: data from the column
        """
        return self.parent._df.iloc[self.index][property]

    def as_rd_mol(self) -> Chem.Mol:
        """
        Get the rdkit molecule object for this instance.
        """
        return Chem.MolFromSmiles(self.smiles)

    def to_file(self, directory, extension=".csv") -> str:
        """
        Write a minimal file containing the SMILES and the ID of the molecule.
        Used for ligrep (.csv is the preferred format).
        """
        filename = os.path.join(directory, self._id + extension)
        if not os.path.isfile(filename):
            with open(filename, "w") as f:
                f.write("SMILES,id\n")
                f.write(f"{self._smiles},{self._id}\n")
        return filename

    def compare(self) -> None:
        """
        Compare the SMILES of the molecule with the original SMILES.
        Visualises side-by-side.
        """
        orig_smiles = self.metadata.original_smiles
        orig_mol = Chem.MolFromSmiles(orig_smiles)
        standardized = Chem.MolFromSmiles(self.smiles)
        print(
            f"Mol ID: {self.id}, Pre/post standardization (structures differ: {Chem.CanonSmiles(orig_smiles) != Chem.CanonSmiles(self.smiles)})"
        )
        return view_difference(orig_mol, standardized)

    def sdf(self) -> List[str] or None:
        """
        Get the SDF file for this molecule.
        """
        sdfs = self.parent._sdf[self.parent._sdf.id == self.id].sdf.values
        if len(sdfs) == 0:
            return None
        else:
            return [decode(sdf) for sdf in sdfs]


class TabularStorage(ChemStore, SMARTSSearchable):
    """
    A simple interface for a chemstore that stores molecules in a pandas dataframe.
    Default columns are: [id, SMILES, original_smiles, metadata, sdf_count, ligprep_metadata]

    The `metadata` column is a dictionary containing standardizer params.
    The `ligprep_metadata` column is a dictionary containing metadata from the ligprep process.

    The store follows the following folder structure:
    """

    def __init__(self, path: str, standardizer=None, identifier=None) -> None:
        super().__init__()

        self._standardizer = (
            PapyrusStandardizer() if standardizer is None else standardizer
        )
        self._identifier = InchiIdentifier() if identifier is None else identifier
        self.path = path = os.path.abspath(path)
        columns = [
            "id",
            "SMILES",
            "original_smiles",
            "library",
            "sdf_count",
            "metadata",
            "~metadata_other",
        ]
        self._df = pd.DataFrame(columns=columns)
        self._mols = dict()
        columns = [
            "library",
            "num_ligands",
            "num_standardized",
            "num_prepared",
            "standardizer_metadata",
            "ligprep_metadata",
        ]
        self._metadata = pd.DataFrame(columns=columns)
        self._sdf = pd.DataFrame(columns=["id", "sdf"])
        self._poses = pd.DataFrame(columns=["id"])
        self._poses.set_index("id", inplace=True)
        self._targets = pd.DataFrame(columns=["id", "pdb"])
        self.load(self.path)

    def summary(self):
        libs = "\n".join([f"\t\t * {lib}" for lib in self._metadata.library.unique()])
        summary = f"""
        Store summary:
        \t * Path: {os.path.abspath(self.path)}
        \t * Number of molecules: {self._df.shape[0]}
        \t * Number of targets: {len(self.targets)}
        \t * Number of molecules with poses: {len(self._poses)}
        \t * Number of libraries: {self._metadata.shape[0]}:
        {libs}
        \t * Standardizer: {self._standardizer.__class__.__name__}
        \t * Identifier: {self._identifier.__class__.__name__}
        \t * SDF files available: {self._sdf.shape[0]}
        """
        print(summary)

    @property
    def standardizer(self) -> ChemStandardizer:
        return self._standardizer

    @property
    def identifier(self):
        return self._identifier

    def load(self, path: str, extension=".tsv"):
        """
        Loads a store from a specified root path.
        """
        # Check if store is a file or a directory
        self.path = path
        store_path = os.path.join(path, f"store{extension}")
        if os.path.isfile(store_path):
            metadata_path = os.path.join(path, f"metadata{extension}")
            self.path = path
            print(f"Loading {self.__class__.__name__} from {path}...")
            self._df = df = load_dataframe(store_path)
            self._metadata = load_dataframe(metadata_path)
            self._sdf = load_dataframe(os.path.join(path, "sdf.tsv"))
            poses_path = os.path.join(path, "poses.tsv")
            try:
                self._poses = load_dataframe(poses_path)
                self._poses.set_index("id", inplace=True)
            except FileNotFoundError:
                warnings.warn(
                    f'Could not find data frame with docked poses: "{poses_path}". Poses will be empty.'
                )
            targets_path = os.path.join(path, "targets.tsv")
            try:
                self._targets = load_dataframe(targets_path)
            except FileNotFoundError:
                warnings.warn(
                    f'Could not find data frame with targets: "{targets_path}". Targets will be empty.'
                )
            print("\tPre-loading molecules...")

            for i, row in tqdm(enumerate(df.iterrows()), total=df.shape[0], ncols=100):
                mol_id = row[1]["id"]
                smiles = row[1]["SMILES"]
                self._mols[mol_id] = TabularMol(self, mol_id, smiles=smiles, index=i)
            self._df["sdf_count"] = self._df["sdf_count"].fillna(0)
        else:
            print(f"Creating new store at <{path}>... ")
        return self

    def save(self, path=None, save_smi_files=False, force=True, **kwargs) -> None:
        """
        Save the store to a file.
        Default is .tsv.
        """
        # TODO: it would be a good idea to make this method thread safe so that it can be called from multiple threads
        if path is None:
            path = self.path

        path = os.path.abspath(path)
        os.makedirs(path, exist_ok=True)
        if os.path.isfile(os.path.join(path, "store.tsv")) and not force:
            print(f"WARNING: Overwriting existing store in <{path}>.")
            if input('Enter "y" to continue...') != "y":
                return

        self._df.to_csv(os.path.join(path, "store.tsv"), index=False, sep="\t")
        self._metadata.to_csv(os.path.join(path, "metadata.tsv"), index=False, sep="\t")
        self._sdf.to_csv(os.path.join(path, "sdf.tsv"), index=False, sep="\t")
        self._poses.to_csv(os.path.join(path, "poses.tsv"), index=True, sep="\t")
        self._targets.to_csv(os.path.join(path, "targets.tsv"), index=False, sep="\t")
        # Create subdirectories
        for subdir in ["config", "sdf", "smiles"]:
            os.makedirs(os.path.join(path, subdir), exist_ok=True)

        ## Saving individual .smi files (for ligprep)
        if save_smi_files:
            smiles_dir = os.path.join(path, "smiles")
            print(f"\nSaving individual .smi files in <{smiles_dir}>...")
            for mol in self.iter_mols():
                mol.to_file(directory=smiles_dir)

    def clear_sdfs(self):
        self._sdf = pd.DataFrame(columns=["id", "sdf"])
        self._df["sdf_count"] = self._df["sdf_count"].fillna(0)

    def clear_poses(self, target: Protein = None):
        if not target:
            self._poses = pd.DataFrame(columns=["id"])
            self._poses.set_index("id", inplace=True)
        else:
            del_col = self._poses.columns[self._poses.columns == f"{target.name}_poses"]
            if len(del_col) > 0:
                del self._poses[del_col[0]]
                self._poses.dropna(axis=0, how="all", inplace=True)

    def apply_standardizer_to_chunk(self, chunk: pd.DataFrame, smiles_col: str):
        smiles = chunk[smiles_col].values
        output = []
        for i, smi in enumerate(smiles):
            try:
                standardized = self.standardizer(smi)[0]
                if standardized is None:
                    raise ValueError(f"Standardizer {self.standardizer} returned None.")
            except Exception as e:
                logging.error(
                    f"Error ({e}) standardizing SMILES: {smi}. "
                    f"Molecule will not be added."
                )
                standardized = None
            output.append((chunk.index[i], standardized, smi))
        return output

    def add_library(
        self,
        library_file,
        smiles_col="SMILES",
        save=False,
        start=0,
        end=-1,
        parallel=False,
        n_jobs=8,
        chunk_size=10,
    ):
        """
        Reads molecules from a file and adds standardized SMILES to the store.

        :param path: path to the library file
        :param smiles_col: name of the column containing the SMILES

        :return: `StoredMol` instance of the added molecule
        """
        print(f"Loading library file... \n\t *{library_file}*")

        added = 0
        if parallel:
            library = load_dataframe(library_file, chunksize=chunk_size)
            para_iter = parallel_iter(
                library, self.apply_standardizer_to_chunk, n_jobs, smiles_col=smiles_col
            )
            output = []
            for chunk in para_iter:
                output.extend(chunk)
            ligand_count = len(output)
            output = [x for x in output if x[1] is not None]
            output = sorted(output, key=lambda x: x[0])
            added = self._add_standardized(
                [x[1] for x in output], [x[2] for x in output], library=library_file
            )
        else:
            library = load_dataframe(library_file)
            ligand_count = len(library)
            try:
                for i, row in tqdm(
                    library.iloc[start:end].iterrows(), total=ligand_count, unit="mol"
                ):
                    smiles = row[smiles_col]
                    mol = self.add_mol_from_smiles(smiles, library=library_file)
                    added += 1 if mol is not None else 0
            except KeyboardInterrupt:
                print("Aborted by user.")

        # Update metadata
        if library_file in self._metadata.library.values:
            self._metadata.loc[
                self._metadata.library == library_file, "num_standardized"
            ] = added
        else:
            new_lib_entry = pd.DataFrame(
                [
                    {
                        "library": library_file,
                        "num_ligands": ligand_count,
                        "num_standardized": added,
                        "num_prepared": 0,
                        "standardizer_metadata": self._standardizer.get_id(),
                        "ligprep_metadata": {},
                    }
                ]
            )
            self._metadata = pd.concat(
                [self._metadata, new_lib_entry], ignore_index=True
            )
        print(
            f"Added {added} ligands to the store ({ligand_count-added} failed/exist already)."
        )
        self.summary()
        if save:
            self.save(save_smi_files=save)

    def _add_standardized(self, standardized: list, original: list, library=None):
        """
        Add a list of standardized SMILES to the store.
        To be used with parallelization.

        :param standardized: list of standardized SMILES
        :param original: list of original SMILES
        :param library: name of the library the molecules belong to

        :return: number of molecules added successfully
        """
        mol_ids, novel, added = [], [], 0
        for i, smiles in tqdm(enumerate(standardized)):
            mol_id = self._identifier(smiles)
            if self._mols.get(mol_id) is None and mol_id not in mol_ids:
                novel += [i]
                mol_ids += [mol_id]
            else:
                logging.warning(
                    f"Molecule already exists in the store or the input library. "
                    f"Duplicate instance ignored: "
                    f"\n\toriginal SMILES: {smiles} "
                    f"\n\tstandardized SMILES: {standardized[i]} "
                    f"\n\tMolecule ID: {mol_id}"
                )

        standardized = [standardized[i] for i in novel]
        original = [original[i] for i in novel]
        numAdded = len(novel)

        library = [library] * numAdded

        data = zip(mol_ids, standardized, original, library)
        new = pd.DataFrame(data, columns=["id", "SMILES", "original_smiles", "library"])
        prev_length = len(self._df)

        self._df = pd.concat([self._df, new], ignore_index=True, axis=0)
        for i, mol_id in enumerate(mol_ids):
            mol = TabularMol(
                self, mol_id, smiles=standardized[i], index=i + prev_length
            )
            self._mols[mol_id] = mol
        return numAdded

    def add_mol(
        self,
        smiles,
        mol_id,
        metadata: dict = None,
        library=None,
        raise_on_existing=False,
        update_existing: Literal["append", "update", None] = None,
    ):
        """
        Add a molecule to the store using its raw SMILES. The SMILES will be standardized and an identifier will be
        calculated.

        :param smiles: SMILES of the molecule to add
        :param mol_id: identifier of the molecule to add
        :param metadata: additional metadata to store with the molecule
        :ligprep_metadata: metadata from the ligprep process
        :param sdf_path: path to the SDF file containing the molecule processed via ligprep
        :param library: name of the library the molecule belongs to
        :param raise_on_existing: whether to raise an error if the molecule already exists in the store
        :param update_existing: whether to update the existing molecule if it already exists in the store

        :return: `StoredMol` instance of the added molecule

        :raises ValueError: if the molecule already exists in the store
        """
        ## Check if mol exists
        if mol_id in self._mols:
            original = self.get_mol(mol_id)
            metadata_orig = original.metadata
            msg = (
                f"Molecule already exists in the store: {mol_id}"
                f"\n\tSMILES={smiles}"
                f"\n\tMetadata_new={metadata}"
                f"\n\tMetadata_existing={metadata_orig}"
            )
            if raise_on_existing:
                raise ValueError(msg)
            if not update_existing:
                msg += "\nMolecule not added, original instance retained. To overwrite, set `update_existing=True`"
                warnings.warn(msg)
                return
            else:
                if update_existing == "update":
                    original.metadata = metadata
                elif update_existing == "append":
                    meta_other = json.loads(metadata_orig["~metadata_other"])
                    if not meta_other:
                        meta_other = []
                    meta_other.append(metadata)
                    metadata_orig["~metadata_other"] = json.dumps(meta_other)
                    original.metadata = metadata_orig
                else:
                    raise ValueError(
                        f"Invalid value for `update_existing`: {update_existing}"
                    )
                return

        ## Add mol (by concatenating to the dataframe)
        original = metadata.pop("original_smiles", None)
        new = pd.DataFrame(
            [
                {
                    "id": mol_id,
                    "SMILES": smiles,
                    "library": library,
                    "original_smiles": original,
                    "metadata": json.dumps(metadata),
                    "~metadata_other": json.dumps([]),
                }
            ]
        )
        self._df = pd.concat([self._df, new], ignore_index=True)
        mol = TabularMol(self, mol_id, smiles=smiles)
        self._mols[mol_id] = mol
        return mol

    def get_mol(self, mol_id) -> TabularMol:
        return self._mols[mol_id]

    def remove_mol(self, mol_id):
        """
        Remove a molecule from the store.
        """
        self._mols.pop(mol_id)
        self._df = self._df[self._df.id != mol_id]

    def add_data(self, mol_id, identifier, data):
        """
        Adds data to a molecule.
        """
        mol = self.get_mol(mol_id)
        mol.add_data(identifier, data)

    def add_sdf(self, sdf_path):
        """
        Adds an SDF file to a molecule.
        """
        sdfs, mol_ids = [], []
        for raw_sdf in sdf_iterator(sdf_path):
            properties = parse_sdf(raw_sdf)
            if len(properties) == 0:
                continue
            mol_id = properties["s_sm_id"]
            mol_ids += [mol_id]
            sdf_mol = "".join(raw_sdf)
            sdfs += [encode(sdf_mol)]
            mol = self.get_mol(mol_id)
            current_count = mol.get_data("sdf_count")
            current_count = 0 if current_count is None else current_count
            mol.add_data("sdf_count", current_count + 1)
        self._sdf = pd.concat(
            [self._sdf, pd.DataFrame({"id": mol_ids, "sdf": sdfs})], ignore_index=True
        )

    def get_mol_ids(self):
        """
        Returns a set of all molecule IDs in the store.
        Good for checking possible overlaps between stores.
        """
        return set(self._mols.keys())

    def get_mol_count(self) -> int:
        return len(self._mols)

    def iter_mols(self) -> Iterator[TabularMol]:
        for mol in self._mols.values():
            yield mol

    def merge_stores(self, other_store):
        """
        Merge another store into this one.

        Not really working yet.
        """
        for mol in other_store.iter_mols():
            self.add_mol(mol.smiles, mol.id, mol.metadata)
        self._metadata = pd.concat(
            [self._metadata, other_store.metadata], ignore_index=True
        )

    def contains_building_block(self, mol_id):
        """
        Returns the set of molecules containing a particular building block.
        """
        pass

    def head(self, n=5):
        return self._df.head(n)

    def __contains__(self, mol_id):
        return self._mols.get(mol_id, None) is not None

    def add_target(self, target: Protein):
        if target.name not in set(self._targets["id"]):
            self._targets = pd.concat(
                [
                    self._targets,
                    pd.DataFrame(
                        {"id": [target.name], "pdb": [encode(target.get_pdb())]}
                    ),
                ],
                ignore_index=True,
            )

    @property
    def targets(self):
        ret = []
        for idx, row in self._targets.iterrows():
            ret.append(
                ProteinFromPDBString(name=row["id"], pdb_string=decode(row["pdb"]))
            )
        return ret

    def add_poses(self, mol_id: str, repr_index: int, poses: Chem.Mol, target: Protein):
        self.add_target(target)
        poses_col = f"{target.name}_poses"
        existing = ""
        if poses_col not in self._poses.columns:
            # Use dtype=object to store compressed hex strings
            self._poses[poses_col] = pd.Series([None] * self._poses.shape[0], dtype=object)
        elif mol_id in self._poses.index:
            previous = self._poses.loc[mol_id, poses_col]
            if previous is not None:
                existing = decode(previous)
        f = io.StringIO()
        poses_sdf = Chem.SDWriter(f)
        repr_smiles = Chem.MolToSmiles(poses, isomericSmiles=True)
        for conf_id in range(poses.GetNumConformers()):
            conf = poses.GetConformer(conf_id)
            name = f"{mol_id}_repr_{repr_index}_target_{target.name}_pose_{conf_id}"
            conf.SetProp("ID", name)
            conf.SetProp("_Name", name)
            conf.SetProp("SMILES_isomeric", repr_smiles)
            props = conf.GetPropsAsDict()
            for prop in props:
                poses.SetProp(prop, str(props[prop]))
            poses_sdf.write(poses, confId=conf_id)
        poses_sdf.flush()
        poses_sdf.close()
        f.flush()
        sdf_str = existing + f.getvalue()
        self._poses.loc[mol_id, poses_col] = encode(sdf_str)
        f.close()

    def get_poses(self, mol_id: str, target: Protein, raw: bool = False):
        poses_col = f"{target.name}_poses"
        if poses_col not in self._poses.columns:
            raise ValueError(f"Poses for target {target.name} not found in store.")
        if mol_id not in self._poses.index:
            raise ValueError(f"No poses for molecule {mol_id} could be found in store.")

        sdf_str = decode(self._poses.loc[mol_id, poses_col])
        if raw:
            return sdf_str
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_str)
        return [mol for mol in suppl]

    def get_complex_for_pose(self, mol_id: str, target: Protein, pose_id: int = 0):
        pose = self.get_poses(mol_id, target)[pose_id]
        # temporary file for pdb string
        pdb_contents = target.get_pdb()
        pdb_file = os.path.join(self.path, "tmp.pdb")
        with open(pdb_file, "w") as f:
            f.write(pdb_contents)
        protein = read_pdb(pdb_file)
        return Chem.CombineMols(protein, pose)

    @staticmethod
    def _find_mols_with_smarts(mols: list[StoredMol], smarts: str):
        return [x for x in mols if smarts_match(x, smarts)]

    def smarts_search(self, smarts: str, n_cpus=None, chunk_size=100):
        n_cpus = n_cpus if n_cpus is not None else os.cpu_count()
        results = []
        for result in parallel_iter(
            self.iter_chunks(chunk_size), self._find_mols_with_smarts, n_cpus, smarts
        ):
            results.extend(result)
        return results

    def get_target(self, name):
        return [x for x in self.targets if x.name == name][0]
