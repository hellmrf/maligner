from rdkit import Chem
from rdkit.Chem import rdFMCS
import datamol as dm

from maligner.mtypes import Mol, MolData


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


def align_and_sanitize_mols(mols: list[Mol], template: Mol) -> list[Mol]:
    """
    Aligns and sanitizes a list of molecules to a given template molecule.

    Args:
        mols (list[Mol]): A list of molecules to be aligned and sanitized.
        template (Mol): The template molecule to align the other molecules to.

    Returns:
        list[Mol]: A list of aligned and sanitized molecules

    """
    template = dm.sanitize_mol(dm.fix_mol(template))
    mols = [dm.align.template_align(mol, template=template) for mol in mols]
    mols = [dm.fix_mol(mol) for mol in mols]
    mols = [dm.sanitize_mol(mol) for mol in mols]
    return mols


def align_moldatas(moldatas: list[MolData], template: Mol) -> list[MolData]:
    """
    Aligns a list of MolData objects to a template molecule.

    Args:
        moldatas (list[MolData]): A list of MolData objects to be aligned.
        template (Mol): The template molecule to align the MolData objects to.

    Returns:
        list[MolData]: The list of aligned MolData objects.
    """
    for moldata in moldatas:
        moldata.mol = dm.align.template_align(mol=moldata.mol, template=template)
        moldata.mol = dm.fix_mol(moldata.mol)
        moldata.mol = dm.sanitize_mol(moldata.mol)

    return moldatas
