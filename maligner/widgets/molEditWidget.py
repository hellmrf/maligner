import sys
from typing import Optional

import numpy as np
from PySide6 import QtCore, QtWidgets
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry.rdGeometry import Point2D

from maligner.mtypes import Mol
from maligner.widgets.molViewWidget import MolWidget

# from rdkit.Chem.AllChem import GenerateDepictionMatching3DStructure


class MolEditWidget(MolWidget):
    def __init__(
        self, mol: Mol, selected_atoms: Optional[list[int]] = None, parent=None
    ):
        if selected_atoms is None:
            selected_atoms = []

        super(MolEditWidget, self).__init__(parent)
        # This sets the window to delete itself when its closed, so it doesn't keep querying the model
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose)

        # Properties
        self._prevmol = None  # For undo
        self.coordlist = None  # SVG coords of the current mols atoms

        self.bondtypes = (
            Chem.rdchem.BondType.names
        )  # A dictionary with all available rdkit bondtypes

        # Points to calculate the SVG to coord scaling
        self.points = [Point2D(0, 0), Point2D(1, 1)]

        # Bind signals to slots

        # When drawing finished, update coordlist of SVG atoms.
        self.finishedDrawing.connect(self.update_coordlist)

        # Init with a mol if passed at construction
        self.mol = mol
        self.selectedAtoms = selected_atoms

    def SVG_to_coord(self, x_svg, y_svg):
        """Function to translate from SVG coords to atom coords using scaling calculated from atomcoords (0,0) and (1,1). Returns rdkit Point2D."""
        if self.drawer is not None:
            scale0 = self.drawer.GetDrawCoords(self.points[0])
            scale1 = self.drawer.GetDrawCoords(self.points[1])

            ax = scale1.x - scale0.x
            bx = scale0.x

            ay = scale1.y - scale0.y
            by = scale0.y

            return Point2D((x_svg - bx) / ax, (y_svg - by) / ay)
        else:
            return Point2D(0.0, 0.0)

    def update_coordlist(self):
        if self.mol is not None:
            self.coordlist = np.array(
                [
                    list(self.drawer.GetDrawCoords(i))
                    for i in range(self.mol.GetNumAtoms())
                ]
            )
            self.logger.debug("Current coordlist:\n%s" % self.coordlist)
        else:
            self.coordlist = None

    def get_nearest_atom(self, x_svg, y_svg):
        if self.mol is not None and self.mol.GetNumAtoms() > 0:
            atomsvgcoords = np.array([x_svg, y_svg])
            # find distance, https://codereview.stackexchange.com/questions/28207/finding-the-closest-point-to-a-list-of-points
            deltas = self.coordlist - atomsvgcoords
            dist_2 = np.einsum("ij,ij->i", deltas, deltas)
            min_idx = np.argmin(dist_2)
            return min_idx, dist_2[min_idx] ** 0.5
        else:
            return None, 1e10  # Return ridicilous long distance so that its not chosen

    def get_nearest_bond(self, x_svg, y_svg):
        if self.mol is not None and self.mol.GetNumAtoms() > 2:
            bondlist = []
            for bond in self.mol.GetBonds():
                bi = bond.GetBeginAtomIdx()
                ei = bond.GetEndAtomIdx()
                avgcoords = np.mean(self.coordlist[[bi, ei]], axis=0)
                bondlist.append(avgcoords)

            bondlist = np.array(bondlist)

            atomsvgcoords = np.array([x_svg, y_svg])
            deltas = bondlist - atomsvgcoords
            dist_2 = np.einsum("ij,ij->i", deltas, deltas)
            min_idx = np.argmin(dist_2)
            return min_idx, dist_2[min_idx] ** 0.5
        else:
            return None, 1e10  # Return ridicilous long distance so that its not chosen

    def get_molobject(self, event):
        """Function that translates coodinates from an event into a molobject"""

        # Recalculate to SVG coords
        viewbox = self.renderer().viewBox()
        size = self.size()
        x = event.pos().x()
        y = event.pos().y()

        # Rescale, divide by the size of the widget, multiply by the size of the viewbox + offset.
        x_svg = float(x) / size.width() * viewbox.width() + viewbox.left()
        y_svg = float(y) / size.height() * viewbox.height() + viewbox.top()
        self.logger.debug("SVG_coords:\t%s\t%s" % (x_svg, y_svg))

        # Identify Nearest atomindex
        atom_idx, atom_dist = self.get_nearest_atom(x_svg, y_svg)
        bond_idx, bond_dist = self.get_nearest_bond(x_svg, y_svg)
        self.logger.debug(
            "Distances to atom %0.2F, bond %0.2F" % (atom_dist, bond_dist)
        )

        # If not below a given threshold, then it was not clicked
        if min([atom_dist, bond_dist]) < 14.0:
            if atom_dist < bond_dist:
                return self.mol.GetAtomWithIdx(int(atom_idx))
            else:
                return self.mol.GetBondWithIdx(int(bond_idx))
        else:
            # Translate SVG to Coords
            return self.SVG_to_coord(x_svg, y_svg)

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            clicked = self.get_molobject(event)
            if isinstance(clicked, Chem.rdchem.Atom):
                self.atom_click(clicked)
            elif isinstance(clicked, Chem.rdchem.Bond):
                self.bond_click(clicked)
            elif isinstance(clicked, Point2D):
                self.logger.debug("Canvas Click")
                self.canvas_click(clicked)

    def atom_click(self, atom):
        """Lookup tables to relate actions to context type with action type"""
        self.select_atom_add(atom)

    def bond_click(self, bond):
        self.select_bond(bond)

    def canvas_click(self, point):
        if len(self.selectedAtoms) > 0:
            self.clearAtomSelection()

    def select_atom(self, atom):
        self.selectAtom(atom.GetIdx())
        # TODO make an unselect atom function

    def select_atom_add(self, atom):
        selidx = atom.GetIdx()
        if selidx in self._selectedAtoms:
            self.unselectAtom(selidx)
        else:
            self.selectAtomAdd(selidx)

    def select_bond(self, bond):
        self.logger.debug("Select_bond not implemented")  # TODO

    def undo(self):
        self.mol = self._prevmol

    def backupMol(self):
        self._prevmol = Chem.Mol(self.mol.ToBinary())


if __name__ == "__main__":
    mol = Chem.MolFromSmiles("CCN(C)C1CCCCC1S")
    rdDepictor.Compute2DCoords(mol)
    myApp = QtWidgets.QApplication(sys.argv)
    molblockview = MolWidget(mol)
    molblockview.show()
    myApp.exec_()
