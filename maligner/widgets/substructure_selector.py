from pathlib import Path
import sys, os

from PySide6 import QtGui, QtWidgets

from rdkit import Chem

from maligner.widgets.molEditWidget import MolEditWidget
from maligner.icons import icon
from maligner.mtypes import MolData


class SubstructureSelectorDialog(QtWidgets.QDialog):

    def __init__(self, moldata: MolData):
        super(SubstructureSelectorDialog, self).__init__()
        self.setWindowTitle(f"Substructure Selector ({moldata.name})")
        self.setWindowIcon(icon('icons8-Cursor.png'))
        self.setModal(True)

        self.editor = MolEditWidget(mol=moldata.mol, selected_atoms=moldata.selected)

        QBtn = QtWidgets.QDialogButtonBox.Cancel | QtWidgets.QDialogButtonBox.Ok

        self.buttonBox = QtWidgets.QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.editor, 0, 0, 1, 1,
                         QtGui.Qt.AlignmentFlag.AlignVCenter | QtGui.Qt.AlignmentFlag.AlignHCenter)
        layout.addWidget(self.buttonBox, 1, 0, 1, 1, QtGui.Qt.AlignmentFlag.AlignRight)
        self.setLayout(layout)

        self.show()


if __name__ == '__main__':
    myApp = QtWidgets.QApplication(sys.argv)
    mol = Chem.MolFromMolFile("molecules/any1.mol")
    # mainWindow = SubstructureSelectorWindow(molecule=mol, loglevel=loglevel)
    dialog = SubstructureSelectorDialog(
        MolData(mol, "any1.mol", Path("molecules/any1.mol"), QtGui.QIcon()))
    dialog.show()
    sys.exit(myApp.exec())
