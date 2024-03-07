from rdkit import Chem
from rdkit.Chem import rdFMCS

from maligner.mtypes import Mol


def find_MCS(mols: list[Mol], flexible=False) -> Mol:
    if len(mols) < 2:
        raise ValueError("At least two molecules are required to find the MCS")

    if flexible:
        mcs = rdFMCS.FindMCS(mols,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             atomCompare=rdFMCS.AtomCompare.CompareAny).smartsString
    else:
        mcs = rdFMCS.FindMCS(mols).smartsString
    mcs_mol = Chem.MolFromSmarts(mcs)
    return mcs_mol


def get_MCS_atoms(mols: list[Mol], flexible=False) -> list[list[int]]:
    if len(mols) < 2:
        raise ValueError("At least two molecules are required to find the MCS")

    mcs_mol = find_MCS(mols, flexible)
    return [list(mol.GetSubstructMatch(mcs_mol)) for mol in mols]
