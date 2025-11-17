import re

import meeko
from rdkit import Chem
from rdkit.Chem import AllChem


def process_stdout(stdout):
    """ Processes the stdout of VinaGPU, returns the affinity of each docking orientation. """
    affinities = []
    is_int = re.compile(r'^\s*\d+\s*$')
    for line in stdout.splitlines():
        if bool(is_int.match(line.decode('utf-8')[:4])):
            orientation_id, affinity, dist1, dist2  = line.split()
            affinities += [float(affinity)]
    return affinities


def prep_mol(mol: Chem.Mol, embed=True, add_hs=True, seed=None):
    Chem.SanitizeMol(mol, Chem.SANITIZE_ALL)
    if add_hs:
        lig = Chem.AddHs(mol)
    else:
        lig = mol
    if embed:
        AllChem.EmbedMolecule(lig, randomSeed=seed if seed is not None else seed)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(lig)
    return meeko_prep.write_pdbqt_string()


def parse_poses(poses):
    pmol = meeko.PDBQTMolecule(poses)
    return meeko.RDKitMolCreate.from_pdbqt_mol(pmol)[0]