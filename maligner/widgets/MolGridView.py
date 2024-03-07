# Just copied from https://stackoverflow.com/a/62031429/5160230, which is an image gallery.
# We can just adapt this for our needs.

import os
import tempfile
from pathlib import Path
from typing import List
from PySide6 import QtCore, QtGui, QtWidgets

from rdkit import Chem
from rdkit.Chem import Draw

from maligner.widgets.substructure_selector import SubstructureSelectorDialog
from maligner.mtypes import Mol, MolData
from maligner import aligner

ICON_SIZE = 200


class StyledItemDelegate(QtWidgets.QStyledItemDelegate):

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.text = option.fontMetrics.elidedText(index.data(), QtCore.Qt.ElideRight, ICON_SIZE)


class MolGridViewWidget(QtWidgets.QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.temp_dir = Path(tempfile.gettempdir())
        # self.choose_btn = QtWidgets.QPushButton(self.tr("Choose"), clicked=self.file_chooser)
        # self.choose_btn.setFixedSize(50, 50)
        self.listview = QtWidgets.QListWidget(
            viewMode=QtWidgets.QListView.IconMode,
            iconSize=ICON_SIZE * QtCore.QSize(1, 1),
            movement=QtWidgets.QListView.Static,
            resizeMode=QtWidgets.QListView.Adjust,
        )
        delegate = StyledItemDelegate(self.listview)
        self.listview.setItemDelegate(delegate)
        self.listview.itemDoubleClicked.connect(self.on_mol_double_click)

        self._filenames: List[str] = []
        self.molecules: List[MolData] = []

        grid_layout = QtWidgets.QGridLayout(self)

        # grid_layout.addWidget(self.choose_btn, 0, 0)
        grid_layout.addWidget(self.listview, 0, 0, 1, 1)

        # self.resize(640, 480)

    @property
    def filenames(self):
        return self._filenames

    @filenames.setter
    def filenames(self, filenames: List[str]):
        self._filenames = filenames
        if len(self._filenames) > 0:
            print(f"{self.load_molecules() = }")

    def on_mol_selection_changed(self, moldata: MolData, selected_atoms: list[int]):
        for moldatas in self.molecules:
            if moldatas == moldata:
                moldatas.selected = selected_atoms
                self.update_moldata_icon(moldatas)
                break
        self.populate_listwidget()

    def on_mol_double_click(self, item: QtWidgets.QListWidgetItem):
        # index = self.listview.currentRow()
        moldata = self.moldata_from_item(item)
        dlg = SubstructureSelectorDialog(moldata)
        dlg.editor.selectionChanged.connect(lambda x: self.on_mol_selection_changed(moldata, x))
        dlg.exec()

    def clear(self):
        self.listview.clear()
        self.molecules = []
        self._filenames = []

    def file_chooser(self):
        self.filenames = QtWidgets.QFileDialog.getOpenFileNames(self, self.tr("Choose molecules"),
                                                                "", "Molecules (*mol)")[0]

    def set_anchor(self, index: int):
        moldata = self.molecules[index]
        old_anchor = moldata.anchor

        for m in self.molecules:
            m.anchor = False

        moldata.anchor = not old_anchor

    def load_molecules(self) -> List[Mol]:
        # TODO: we'll need to handle filetype here!
        if not self.filenames:
            return []

        self.molecules = []

        for file in self.filenames:
            name = os.path.basename(file)
            mol = Chem.MolFromMolFile(file)
            if mol is None:
                raise ValueError(f"Could not read mol file {file}")
            qicon = self.molecule_to_icon(mol, name)
            moldata = MolData(mol, name, Path(file), qicon)
            self.molecules.append(moldata)

        self.populate_listwidget()
        return self.molecules

    def update_all_moldata_icons(self):
        for moldata in self.molecules:
            self.update_moldata_icon(moldata)

    def update_moldata_icon(self, moldata: MolData) -> None:
        img_file = self.temp_dir / f"{moldata.name}.png"

        Draw.MolToImageFile(moldata.mol,
                            filename=img_file,
                            size=(ICON_SIZE, ICON_SIZE),
                            highlightAtoms=moldata.selected)
        qimg = QtGui.QImage(img_file)
        pixmap = QtGui.QPixmap.fromImage(qimg)
        moldata.qicon = QtGui.QIcon(pixmap)

    def molecule_to_icon(self, mol: Mol, name: str) -> QtGui.QIcon:
        if not name.endswith(".mol"):
            name += ".mol"

        img_file = self.temp_dir / f"{name}.png"

        Draw.MolToImageFile(mol, filename=img_file, size=(ICON_SIZE, ICON_SIZE))
        qimg = QtGui.QImage(img_file)
        pixmap = QtGui.QPixmap.fromImage(qimg)
        return QtGui.QIcon(pixmap)

    def moldata_from_item(self, item: QtWidgets.QListWidgetItem) -> MolData:
        index = self.listview.indexFromItem(item).row()
        return self.molecules[index]

    def populate_listwidget(self):
        self.listview.clear()

        for moldata in self.molecules:
            name = f"⚓ {moldata.name}" if moldata.anchor else moldata.name
            it = QtWidgets.QListWidgetItem(name)
            it.setIcon(moldata.qicon)
            self.listview.addItem(it)

    def compute_MCS(self):
        selected_atoms = aligner.get_MCS_atoms([m.mol for m in self.molecules])
        for i, moldata in enumerate(self.molecules):
            moldata.selected = selected_atoms[i]
            self.update_moldata_icon(moldata)

    def run_alignment(self):
        anchor_mol = next((m for m in self.molecules if m.anchor), None)
        if anchor_mol is None:
            QtWidgets.QMessageBox.warning(self, "No anchor molecule",
                                          "Please select an anchor molecule for alignment.")
            return

        mcs = aligner.find_MCS([m.mol for m in self.molecules])
        self.molecules = aligner.align_moldatas(self.molecules, mcs)
        self.update_all_moldata_icons()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = MolGridViewWidget()
    w.show()
    sys.exit(app.exec())
