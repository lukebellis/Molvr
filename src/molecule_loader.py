import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem  # Import AllChem for generating coordinates

def load_molecule(mol_file_path):
    """
    Load a molecule from a .mol file and generate its 3D coordinates.
    
    Parameters:
    - mol_file_path: str, path to the .mol file.
    
    Returns:
    - The molecule object if successful, None otherwise.
    """
    mol = Chem.MolFromMolFile(mol_file_path)
    if mol is not None:
        print(f"Molecule loaded successfully from {mol_file_path}.")
        # Generate 3D coordinates for better visualization. For 2D, use Compute2DCoords.
        AllChem.Compute2DCoords(mol)
        return mol
    else:
        print(f"Failed to load molecule from {mol_file_path}.")
        return None

def load_molecules_from_directory(directory_path):
    """
    Load all .mol files from a specified directory and generate their 3D coordinates.
    
    Parameters:
    - directory_path: str, the path to the directory containing .mol files.
    
    Returns:
    - A list of molecule objects loaded from the directory.
    """
    molecules = []
    # Construct the full path for .mol file search
    search_path = os.path.join(directory_path, "*.mol")
    for mol_file in glob.glob(search_path):
        mol = load_molecule(mol_file)
        if mol is not None:
            molecules.append(mol)
    return molecules

def get_atom_positions(mol):
    """
    Returns a list of tuples containing the 3D positions of each atom in the molecule.
    
    Parameters:
    - mol: The RDKit molecule object.
    
    Returns:
    - A list of (x, y, z) tuples for each atom.
    """
    positions = mol.GetConformer().GetPositions()
    return positions

def get_bonds(mol):
    """
    Returns a list of tuples indicating atom indices that form bonds.
    
    Parameters:
    - mol: The RDKit molecule object.
    
    Returns:
    - A list of (start_atom_idx, end_atom_idx) tuples for each bond.
    """
    bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]
    return bonds
