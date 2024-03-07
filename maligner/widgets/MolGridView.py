# Just copied from https://stackoverflow.com/a/62031429/5160230, which is an image gallery.
# We can just adapt this for our needs.

import os
import tempfile
from pathlib import Path
from typing import List, TypeAlias
from PySide6 import QtCore, QtGui, QtWidgets
from typing import Optional

from rdkit import Chem
from rdkit.Chem import Draw
import datamol as dm

Mol: TypeAlias = dm.Mol

ICON_SIZE = 200


class StyledItemDelegate(QtWidgets.QStyledItemDelegate):

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.text = option.fontMetrics.elidedText(index.data(), QtCore.Qt.ElideRight, ICON_SIZE)


class MolGridViewWidget(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.temp_dir = Path(tempfile.gettempdir())
        self.choose_btn = QtWidgets.QPushButton(self.tr("Choose"),
                                                clicked=self.on_choose_btn_clicked)
        self.choose_btn.setFixedSize(50, 50)
        self.path_le = QtWidgets.QLineEdit()
        # self.back_btn = QtWidgets.QPushButton(self.tr("^"), clicked=self.on_back_btn_clicked)
        # self.back_btn.setFixedSize(20, 20)
        self.pixmap_lw = QtWidgets.QListWidget(
            viewMode=QtWidgets.QListView.IconMode,
            iconSize=ICON_SIZE * QtCore.QSize(1, 1),
            movement=QtWidgets.QListView.Static,
            resizeMode=QtWidgets.QListView.Adjust,
        )
        delegate = StyledItemDelegate(self.pixmap_lw)
        self.pixmap_lw.setItemDelegate(delegate)

        self._filenames: List[str] = []
        self.molecules: List[Mol] = []

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        grid_layout = QtWidgets.QGridLayout(central_widget)

        grid_layout.addWidget(self.choose_btn, 0, 0, 2, 1)
        grid_layout.addWidget(self.path_le, 0, 1)
        # grid_layout.addWidget(self.back_btn, 1, 1, alignment=QtCore.Qt.AlignRight)
        grid_layout.addWidget(self.pixmap_lw, 2, 0, 1, 2)

        self.resize(640, 480)

        # self.timer_loading = QtCore.QTimer(interval=50, timeout=self.load_image)
        # self.filenames_iterator = None

    @property
    def filenames(self):
        return self._filenames

    @filenames.setter
    def filenames(self, filenames: List[str]):
        self._filenames = filenames
        if len(self._filenames) > 0:
            print(f"{self.load_molecules() = }")

    @QtCore.Slot()
    def on_choose_btn_clicked(self):
        self.filenames = QtWidgets.QFileDialog.getOpenFileNames(self, self.tr("Choose molecules"),
                                                                "", "Molecules (*mol)")[0]

    def load_molecules(self) -> List[Mol]:
        # TODO: we'll need to handle filetype here!
        if not self.filenames:
            return []
        self.molecules = [Chem.MolFromMolFile(file) for file in self.filenames]
        self.populate_listwidget()
        return self.molecules

    def molecule_to_icon(self, mol: Mol, name: str) -> QtGui.QIcon:
        if not name.endswith(".mol"):
            name += ".mol"

        img_file = self.temp_dir / f"{name}.png"

        Draw.MolToImageFile(mol, filename=img_file, size=(ICON_SIZE, ICON_SIZE))
        qimg = QtGui.QImage(img_file)
        pixmap = QtGui.QPixmap.fromImage(qimg)
        return QtGui.QIcon(pixmap)

    # def molecules_to_icons(self) -> List[QtGui.QIcon]:
    #     icons = [self.molecule_to_icon(mol) for mol in self.molecules]
    #     return icons

    def render_mol_previews(self):
        pass

    def populate_listwidget(self):
        for i, mol in enumerate(self.molecules):
            name = os.path.basename(self.filenames[i])
            icon = self.molecule_to_icon(mol, name)
            # smiles = dm.to_smiles(mol)
            it = QtWidgets.QListWidgetItem(name)
            it.setIcon(icon)
            self.pixmap_lw.addItem(it)


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = MolGridViewWidget()
    w.show()
    sys.exit(app.exec())
