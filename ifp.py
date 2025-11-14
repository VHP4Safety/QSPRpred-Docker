import json
import os
from random import choice
from string import ascii_uppercase

import numpy as np
import pandas as pd
import prolif as plf
from qsprpred.data import MoleculeTable
from qsprpred.data.descriptors.sets import DataFrameDescriptorSet
from rdkit import Chem, DataStructs
from spock.docking.vina.cpu_local import VinaDockingCPULocal, VinaProtein
from spock.prep import Dimorphite, Ligprep
from spock.storage import TabularStorage
from spock.utils.standardizers.papyrus import PapyrusStandardizer


def initialize(lib_name, lib_path):
    # load molecules
    standardizer = PapyrusStandardizer()
    standardizer.small_molecule_min_mw = 50
    store = TabularStorage(path=lib_name, standardizer=standardizer)
    store.add_library(lib_path, parallel=True, chunk_size=10, n_jobs=5, save=True)
    return store

def load(lib_name):
    # load molecules
    standardizer = PapyrusStandardizer()
    standardizer.small_molecule_min_mw = 50
    store = TabularStorage(path=lib_name, standardizer=standardizer)
    return store


def add_Hs(store, prep):
    # add hydrogens
    if prep == "Schrodinger":
        ligprep = Ligprep(schrodinger_path=None)
        ligprep.params.num_stereoisomers = 4
        ligprep.process_store(store, workers=4, wait=True, chunks=2)
    else:
        ligprep = Dimorphite()
        ligprep.max_variants = 1
        ligprep.process_store(store)

    store.save(force=True)

    return store


def dock(store, protein_id):
    # dock
    protein = VinaProtein(
        f"{protein_id}",  # just a name
        f"./models/{protein_id}_final/{protein_id}.pdb",
        f"./models/{protein_id}_final/{protein_id}.pdbqt"
    )

    center_coords = {"8t6j_b": [204, 170, 253] #[194, 170, 255]
                     , "4or2": [165, 178, 189]}

    docking = VinaDockingCPULocal(
        protein=protein,
        n_cpus=12,
        box_spec={
            "center": center_coords[protein_id],
            "box_size": [26, 26, 26]
        },
        embed_mols=True,  # set to False if conformers are already generated
        exhaustiveness=8,
        seed=42,
    )
    docking.dock_storage(store, chunk_size=1, overwrite=True, save=True)
    return store


def get_IFPs(store, protein_id):
    # get IFPs
    protein_file = str(f"./models/{protein_id}_final/{protein_id}.pdb")
    rdkit_prot = Chem.MolFromPDBFile(protein_file, removeHs=False)
    protein_mol = plf.Molecule(rdkit_prot)

    fp = plf.Fingerprint(
        [
            "Anionic", "CationPi", "Cationic", "EdgeToFace", "FaceToFace", "HBAcceptor",
            "HBDonor", "Hydrophobic", "MetalAcceptor", "MetalDonor", "PiCation",
            "PiStacking", "VdWContact", "XBAcceptor", "XBDonor"
        ],
        count=True
    )

    df = store._df.iloc[range(len(store._df))]
    mol_id_lookup = pd.Series(df.original_smiles.values, index=df["id"]).to_dict()

    df_comb = pd.DataFrame()
    for i, mol_id in enumerate(store.get_mol_ids()):
        try:
            poses = store.get_poses(mol_id, target=store.targets[0])
            with Chem.SDWriter(f"{mol_id_lookup[mol_id]}_{protein_id}.sdf") as w:
                for pose in poses:
                    w.write(pose)
        except ValueError:
            continue
        scores = [x.GetProp("vina_energy_total") for x in poses]
        pose_iterable = [plf.Molecule(x) for x in poses]

        fp.run_from_iterable(pose_iterable, protein_mol)
        df = fp.to_dataframe(index_col="Pose", drop_empty=False)
        df = pd.DataFrame(df.values, columns=[f"ifp_{x[1]}_{x[2]}" for x in df.columns])
        df["mol_ID"] = mol_id
        df["pose_ID"] = [mol_id + f"_{i}" for i in range(len(poses))]
        df["vina_energy_total"] = scores
        df_comb = pd.concat([df_comb, df], axis=0)
    # Convert df to bitvectors
    ids = df_comb["mol_ID"].to_list()
    scores = df_comb["vina_energy_total"].astype(float).to_list()

    # Select top three poses per compound
    df = pd.DataFrame(
        {
            "mol_ID": ids,
            "pose_ID": df_comb["pose_ID"].to_list(),
            "vina_energy_total": df_comb["vina_energy_total"].astype(float).to_list(),
        }
    )
    df = df.groupby("mol_ID", group_keys=False
                   ).apply(lambda x: x.nsmallest(3, "vina_energy_total"))
    df_comb = df_comb[df_comb["pose_ID"].isin(df["pose_ID"])]

    # Average fingerprint bits for top three poses
    df_comb = df_comb.drop(columns=["pose_ID", "vina_energy_total"])
    df_comb[df_comb.columns.difference(["mol_ID"])
           ] = df_comb.groupby("mol_ID")[df_comb.columns.difference(["mol_ID"]
                                                                   )].transform("mean")

    # Keep one row for each unique compound
    df_comb_profile = df_comb.drop_duplicates(subset="mol_ID",
                                              keep="first").drop("mol_ID", axis=1)

    # turn into bitvectors
    bitvectors = plf.to_bitvectors(df_comb_profile)

    # Add fingerprint bits to numpy array
    num_bits = bitvectors[0].GetNumBits()
    arr = np.zeros((len(bitvectors), num_bits), dtype=np.int8)

    for i, bitvector in enumerate(bitvectors):
        DataStructs.ConvertToNumpyArray(bitvector, arr[i])

    bits = pd.DataFrame(arr, columns=df_comb_profile.columns)

    with open("./models/interactions.json") as json_file:
        interactions = json.load(json_file)
    bits_aligned = bits.reindex(columns=interactions[protein_id], fill_value=0)
    bits_aligned["ID"] = [x for x in list(store.get_mol_ids()) if x in ids]

    return bits_aligned


def predict(smiles_list, model, protein_id, prep):
    if not os.path.exists('./poses'):
        os.makedirs('./poses')

    mol_id = ''.join(choice(ascii_uppercase) for i in range(12))

    lib_name = f"./poses/{protein_id + '_' + mol_id}"

    df = pd.DataFrame({"SMILES": smiles_list})
    lib_path = f"./poses/{mol_id}.csv"

    df.to_csv(lib_path, index=False)

    store = initialize(lib_name, lib_path)
    df = store._df.iloc[range(len(store._df))]
    original_smiles = pd.Series(df["id"].values,index=df.original_smiles).to_dict()
    original_order = [original_smiles[smile] for smile in smiles_list]
    store = add_Hs(store, prep)

    store = dock(store, protein_id)

    bits_aligned = get_IFPs(store, protein_id)
    bits_aligned = bits_aligned.drop(columns=["ID"])


    df = store._df.iloc[range(len(smiles_list))]
    file_name = "table"
    mt = MoleculeTable(
        df=df, store_dir="./poses", name=f"{file_name}", overwrite=True, random_state=42
    )
    ifp = bits_aligned.set_index(mt.getDF().index)

    # Add descriptors to dataset
    descriptors = DataFrameDescriptorSet(ifp)
    mt.addDescriptors([descriptors], recalculate=True)
    new_order = mt.getDF()["id"].tolist()
    res = [new_order.index(idx) for idx in original_order]


    predictions = model.predictProba(mt.getDescriptors().to_numpy())

    predictions[0] = predictions[0][res]
    predictions = [[1] if x > 0.5 else [0] for x in predictions[0][:,1]]

    series_index = mt.getDF().index.tolist()
    series_index_res = [series_index[x] for x in res]

    within_ad = model.applicabilityDomain.contains(mt.getDescriptors()).reindex(series_index_res)

    return predictions, within_ad, lib_name