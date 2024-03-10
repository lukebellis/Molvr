import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
from molecule_processor import load_molecules_from_directory, get_atom_positions, get_bonds

import numpy as np

def draw_atom(position, element):
    """Draws an atom (as a sphere) at the given position."""
    glColor3fv((1, 0, 0))  # Color, you might want to change colors based on the element
    glPushMatrix()
    glTranslate(*position)
    glutSolidSphere(0.5, 20, 20)  # Sphere for the atom, adjust size as needed
    glPopMatrix()

def draw_bond(start_pos, end_pos):
    """Draws a bond (as a line) between two atoms."""
    glColor3fv((1, 1, 1))  # Color for the bond
    glBegin(GL_LINES)
    glVertex3fv(start_pos)
    glVertex3fv(end_pos)
    glEnd()

def visualize_molecule(molecule_path):
    pygame.init()
    display = (800, 600)
    pygame.display.set_mode(display, DOUBLEBUF|OPENGL)
    gluPerspective(45, (display[0]/display[1]), 0.1, 50.0)
    glTranslatef(0.0, 0.0, -30)  # Adjust the translation to fit the molecule in view

    mol = load_molecule(molecule_path)
    if mol is None:
        return

    atom_positions = get_atom_positions(mol)
    bonds = get_bonds(mol)

    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                quit()

        glRotatef(1, 3, 1, 1)
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

        for position in atom_positions:
            draw_atom(position, 'C')  # Assuming carbon for simplicity, adjust as needed

        for start_idx, end_idx in bonds:
            start_pos, end_pos = atom_positions[start_idx], atom_positions[end_idx]
            draw_bond(start_pos, end_pos)

        pygame.display.flip()
        pygame.time.wait(10)

if __name__ == "__main__":
    molecule_path = 'assets/molecules/your_molecule.mol'  # Adjust path to your molecule
    visualize_molecule(molecule_path)
