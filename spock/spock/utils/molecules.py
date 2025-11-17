from rdkit import Chem

from spock.storage.base import StoredMol


def smarts_match(mol: StoredMol | Chem.Mol, smarts):
    """
    Returns True if the molecule matches the SMARTS pattern, False otherwise.
    """
    if isinstance(mol, StoredMol):
        mol_rd = mol.as_rd_mol()
        patt = Chem.MolFromSmarts(smarts)
        has_match  = mol_rd.HasSubstructMatch(patt)
        if not has_match:
            for mol in mol.representations:
                if mol.HasSubstructMatch(patt):
                    has_match = True
                    break
        return has_match
    else:
        return mol.HasSubstructMatch(patt)
