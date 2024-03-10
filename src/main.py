from src.molecule_loader import load_molecules_from_directory

directory_path = 'assets/molecules/'  # Adjust as necessary

# Load all molecules from the directory
molecules = load_molecules_from_directory(directory_path)

# Example: Process the first molecule (ensure the list is not empty)
if molecules:
    mol = molecules[0]
    atom_positions = get_atom_positions(mol)
    bonds = get_bonds(mol)
    
    # Now, atom_positions and bonds are ready for use in visualization
    print(atom_positions)
    print(bonds)
