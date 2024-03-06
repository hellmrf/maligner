from pathlib import Path
import sys, os

from typing import Optional
from PySide6 import QtGui, QtWidgets
from rdkit import Chem

from maligner.widgets.molEditWidget import MolEditWidget


class SubstructureSelectorWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, filename: Optional[Path | str], loglevel="WARNING"):
        super(SubstructureSelectorWindow, self).__init__()
        self.loglevels = ["Critical", "Error", "Warning", "Info", "Debug", "Notset"]
        self.editor = MolEditWidget()
        self._filename = filename
        self.init_GUI(filename=filename)
        # TODO: conectar signal selectionChanged
        # self.editor.selectionChanged.connect(self.setAtomTypeName)
        self.editor.logger.setLevel(loglevel)

    def get_pixmap(self, name: str) -> str:
        p = Path(__file__).parent.parent / 'pixmaps' / name
        return QtGui.QPixmap(p)

    #Properties
    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        if filename != self._filename:
            self._filename = filename
            self.setWindowTitle(str(filename))

    def init_GUI(self, filename=None):
        self.setWindowTitle("Substructure selector")
        self.setWindowIcon(QtGui.QIcon(self.get_pixmap('icons8-Cursor.png')))
        self.setGeometry(100, 100, 200, 150)

        self.center = self.editor
        self.center.setFixedSize(600, 600)
        self.setCentralWidget(self.center)
        self.filename = filename

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.setup_components()

        self.infobar = QtWidgets.QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        if self.filename is not None:
            self.editor.logger.info("Loading model from %s" % self.filename)
            self.load_mol_file(filename)

        self.editor.sanitizeSignal.connect(self.infobar.setText)
        self.show()

    # Function to setup status bar, central widget, menu bar, tool bar
    def setup_components(self):
        self.myStatusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)

        self.create_actions()
        self.create_tool_bars()

    def create_tool_bars(self):
        self.mainToolBar = self.addToolBar('Main')
        #Main action bar
        self.mainToolBar.addAction(self.cleanCoordinatesAction)
        self.mainToolBar.addAction(self.undoAction)
        self.mainToolBar.addAction(self.deleteMoleculeAction)

    def load_mol_file(self, filename):
        self.filename = filename
        mol = Chem.MolFromMolFile(str(self.filename), sanitize=False, strictParsing=False)
        self.editor.mol = mol
        self.statusBar().showMessage("File opened")

    def openFile(self):
        self.filename, self.filterName = QtWidgets.QFileDialog.getOpenFileName(
            self, caption="Open MOL file", filter=self.filters)
        self.load_mol_file(self.filename)

    def save_file(self):
        if self.filename != None:
            Chem.MolToMolFile(self.editor.mol, str(self.filename))
        else:
            self.save_as_file()

    def save_as_file(self):
        self.filename, self.filterName = QtWidgets.QFileDialog.getSaveFileName(self,
                                                                               filter=self.filters)
        if (self.filename != ''):
            if self.filename[-4:].upper() != ".MOL":
                self.filename = self.filename + ".mol"
            Chem.MolToMolFile(self.editor.mol, str(self.filename))
            # file = open(self.fileName, 'w')
            # file.write(self.textEdit.toPlainText())
            self.statusBar().showMessage("File saved", 2000)

    def clear_canvas(self):
        self.editor.mol = None
        self.filename = None
        self.exit_file()
        self.statusBar().showMessage("Molecule removed from canvas", 2000)

    def closeEvent(self, event):
        self.exit_file()
        event.ignore()

    def exit_file(self):
        exit(0)  #TODO, how to exit QtWidgets.QApplication from within class instance?

    def about_help(self):
        QtWidgets.QMessageBox.about(
            self, "About Simple Molecule Editor",
            """A Simple Molecule Editor where you can edit molecules\nBased on RDKit! http://www.rdkit.org/ \nSome icons from http://icons8.com\n\nSource code: https://github.com/EBjerrum/rdeditor"""
        )

    def set_log_level(self):
        loglevel = self.sender().objectName().split(':')[-1].upper()
        self.editor.logger.setLevel(loglevel)

    # Function to create actions for menus and toolbars
    def create_actions(self):
        self.exitAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap('icons8-Shutdown.png')),
                                        'E&xit',
                                        self,
                                        shortcut="Ctrl+Q",
                                        statusTip="Exit the Application",
                                        triggered=self.exit_file)

        #Misc Actions
        self.undoAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap('prev.png')),
                                        'U&ndo',
                                        self,
                                        shortcut="Ctrl+Z",
                                        statusTip="Undo/Redo changes to molecule Ctrl+Z",
                                        triggered=self.editor.undo,
                                        objectName="undo")

        self.deleteMoleculeAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap('icons8-Trash.png')),
                                                  'Delete molecule from set (Del)',
                                                  self,
                                                  shortcut="Del",
                                                  statusTip="Remove this molecule from set (Del)",
                                                  triggered=self.clear_canvas,
                                                  objectName="Clear Canvas")

        self.cleanCoordinatesAction = QtGui.QAction(
            QtGui.QIcon(self.get_pixmap('icons8-Broom.png')),
            'Recalculate coordinates (F9)',
            self,
            shortcut="F9",
            statusTip="Re-calculates coordinates and redraw (F9)",
            triggered=self.editor.canon_coords_and_draw,
            objectName="Recalculate Coordinates")


class SubstructureSelectorDialog(QtWidgets.QDialog):

    def __init__(self, filename: Optional[Path | str], loglevel="WARNING"):
        super().__init__()
        self.ss_window = SubstructureSelectorWindow(filename, loglevel)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.ss_window)
        self.setLayout(layout)
        self.setWindowTitle("Substructure Selector")
        self.setWindowIcon(QtGui.QIcon(self.ss_window.get_pixmap('icons8-Cursor.png')))


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    myApp = QtWidgets.QApplication(sys.argv)
    argv1 = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if argv1 and argv1.is_file():
        mainWindow = SubstructureSelectorWindow(filename=argv1, loglevel=loglevel)
    else:
        mainWindow = SubstructureSelectorWindow(filename=None, loglevel=loglevel)
    sys.exit(myApp.exec())


if __name__ == '__main__':
    launch(loglevel="DEBUG")
