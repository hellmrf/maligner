from rdkit import Chem
from rdkit.Chem import rdFMCS

from maligner.mtypes import Mol


def find_MCS(mols: list[Mol], flexible=False) -> Mol:
    """ Finds the Maximum Common Substructure (MCS) among a list of molecules.

    Args:
        mols (list[Mol]): A list of RDKit Mol objects representing the input molecules.
        flexible (bool, optional): Specifies whether to allow flexible matching or not.
            If set to True, matching structures with different atoms or with different 
            bond types will be allowed. For example, a carbon atom will match with a
            nitrogen atom, and a single bond will match with a double bond.
            Defaults to False.

    Returns:
        Mol: The RDKit Mol object representing the MCS.
    """
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
    """
    Returns the atoms that correspond to the Maximum Common Substructure (MCS) in a list of molecules.

    Args:
        mols (list[Mol]): A list of RDKit Mol objects representing the molecules.
        flexible (bool, optional): Whether to allow flexible matching of the MCS. 
            See `maligner.aligner.find_MCS` for more info. Defaults to False.

    Returns:
        list[list[int]]: A list of lists, where each inner list contains the atom indices that correspond to the MCS in each molecule.

    Raises:
        ValueError: If less than two molecules are provided.
    """
    if len(mols) < 2:
        raise ValueError("At least two molecules are required to find the MCS")

    mcs_mol = find_MCS(mols, flexible)
    return [list(mol.GetSubstructMatch(mcs_mol)) for mol in mols]
