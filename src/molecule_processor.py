from rdkit import Chem
from rdkit.Chem import AllChem

def load_molecule(mol_file_path):
    """
    Load a molecule from a .mol file and optionally generate 3D coordinates.
    
    Parameters:
    - mol_file_path: str, path to the .mol file.
    
    Returns:
    - The molecule object if successful, None otherwise.
    """
    mol = Chem.MolFromMolFile(mol_file_path)
    if mol is not None:
        print(f"Molecule loaded successfully from {mol_file_path}.")
        # Generate 2D coordinates for visualization purposes.
        # Use AllChem.Compute2DCoords for 2D or AllChem.Compute3DCoords for 3D
        AllChem.Compute2DCoords(mol)
        return mol
    else:
        print(f"Failed to load molecule from {mol_file_path}.")
        return None

def get_atom_positions(mol):
    """
    Extracts the positions of each atom in a molecule.
    
    Parameters:
    - mol: The RDKit molecule object.
    
    Returns:
    - A list of (x, y, z) tuples for the position of each atom.
    """
    positions = mol.GetConformer().GetPositions()
    return positions

def get_bonds(mol):
    """
    Extracts the bonding information from a molecule.
    
    Parameters:
    - mol: The RDKit molecule object.
    
    Returns:
    - A list of tuples, each containing the indices of two atoms that form a bond.
    """
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    return bonds
