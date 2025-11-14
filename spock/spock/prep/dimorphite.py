import math
import os

import pandas as pd
from dimorphite_dl import protonate_smiles
from rdkit import Chem
from tqdm.auto import tqdm

from spock.prep.base import LigandPreparation


class Dimorphite(LigandPreparation):

    def __init__(self, ph_range=(6, 7), max_variants=128, pka_precision=0.5):
        super().__init__()
        self.ph_range = ph_range
        self.max_variants = max_variants
        self.pka_precision = pka_precision
        self._args = {
            "ph_min":self.ph_range[0],
            "ph_max":self.ph_range[1],
            "max_variants":self.max_variants,
            "label_states":False,
            "precision":self.pka_precision
        }

    def run(self, filename):
        """
        Run ligprep on a single molecule
        """
        filename = os.path.abspath(filename)
        input_extension = filename.split(".")[-1]
        in_file = os.path.basename(filename).split(".")[0]
        output_filename = os.path.join(os.path.dirname(filename), f"{in_file}.sdf")
        formats = ["csv"]
        if input_extension not in formats:
            raise ValueError(f"Unsupported input format: {input_extension}")
        df_in = pd.read_csv(filename)
        sdf_writer = Chem.SDWriter(output_filename)
        for idx, row in df_in.iterrows():
            mol_id = row.id
            smiles = row.SMILES
            protonated = protonate_smiles(smiles, **self._args)
            for smile in protonated:
                mol = Chem.MolFromSmiles(smile)
                mol.SetProp("_Name", mol_id)
                mol.SetProp("s_sm_id", mol_id)
                sdf_writer.write(mol, 0)
        sdf_writer.close()
        return output_filename

    def process_store(self, store, cpus=None, chunks=1):
        """
        Process all molecules in a store
        """
        new_df = store._df.copy()
        new_df.drop(columns=["metadata", "original_smiles", "library", "sdf_count"],
                    inplace=True)
        ## Make 'SMILES' the first column
        cols = new_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        new_df = new_df[cols]
        new_df["NAME"] = new_df["id"]
        # Split the csv file into chunks
        mols_per_chunk = math.ceil(len(new_df) / chunks)
        store.clear_sdfs()
        for i in tqdm(range(chunks), desc="Processing files...", total=chunks):
            chunk = new_df.iloc[i * mols_per_chunk:(i + 1) * mols_per_chunk]
            prefix = "store_ligprep"
            path = os.path.join(store.path, f"{prefix}.csv")
            chunk[["SMILES", "NAME", "id"]].to_csv(path, index=False, header=True)
            output_sd_path = self.run(path)
            store.add_sdf(output_sd_path)
            if os.path.exists(output_sd_path):
                os.remove(output_sd_path)
            if os.path.exists(path):
                os.remove(path)
        store.save(force=True)
