import sys, os
from pathlib import Path
from typing import List, Optional

from rdkit import Chem
from PySide6 import QtGui, QtWidgets, QtCore
from maligner.mtypes import MolData
from maligner.widgets.MolGridView import MolGridViewWidget

from maligner.widgets.molEditWidget import MolEditWidget
from maligner.widgets.substructure_selector import SubstructureSelectorDialog


class MainWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, filenames: Optional[list[Path | str]] = None):
        super(MainWindow, self).__init__()
        self.loglevels = ["Critical", "Error", "Warning", "Info", "Debug", "Notset"]
        # self.editor = MolEditWidget()
        # self.substructure_selector = SubstructureSelectorDialog()
        self._filenames = filenames
        self.molgridview = MolGridViewWidget()
        self.init_GUI()
        # TODO: selectionChanges ainda não existe
        # self.substructure_selector.selectionChanged.connect(self.setAtomTypeName)
        # self.editor.logger.setLevel(loglevel)

    def get_pixmap(self, name: str) -> QtGui.QPixmap:
        p = Path(__file__).parent / 'pixmaps' / name
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

    def init_GUI(self):
        self.setWindowTitle(r"maligner |\ An Open-Source Molecular Alignment Tool")
        self.setWindowIcon(QtGui.QIcon(self.get_pixmap("appicon.svg.png")))
        self.setGeometry(400, 400, 700, 500)

        self.setCentralWidget(self.molgridview)

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.setup_components()

        self.infobar = QtWidgets.QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        self.show()

    # Function to setup status bar, central widget, menu bar, tool bar
    def setup_components(self):
        self.myStatusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)

        self.create_actions()
        self.create_menus()
        self.create_tool_bars()

    # Actual menu bar item creation
    def create_menus(self):
        self.fileMenu = self.menuBar().addMenu("&File")
        self.toolMenu = self.menuBar().addMenu("&Tools")
        self.helpMenu = self.menuBar().addMenu("&Help")

        # File
        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addAction(self.saveAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        # Tools
        self.toolMenu.addAction(self.anchorAction)
        self.toolMenu.addAction(self.deleteMoleculeAction)
        self.toolMenu.addAction(self.openSelectorAction)

        #Help menu
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addAction(self.aboutQtAction)
        #Debug level sub menu
        self.loglevelMenu = self.helpMenu.addMenu("Logging Level")
        for loglevel in self.loglevels:
            self.loglevelMenu.addAction(self.loglevelactions[loglevel])

    def create_tool_bars(self):
        self.mainToolBar = self.addToolBar('Main')
        #Main action bar
        self.mainToolBar.addAction(self.openAction)
        self.mainToolBar.addAction(self.saveAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.anchorAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.deleteMoleculeAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.openSelectorAction)

    def load_mol_file(self, filename):
        self.filename = filename
        mol = Chem.MolFromMolFile(str(self.filename), sanitize=False, strictParsing=False)
        self.statusBar().showMessage("File opened")

    def open_file(self):
        self.molgridview.file_chooser()
        self.molgridview.load_molecules()

    def saveFile(self):
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
            self.statusBar().showMessage("File saved", 2000)

    def open_selector(self):
        citem = self.molgridview.listview.currentItem()
        if citem is None:
            self.statusBar().showMessage("No molecule selected", 2000)
            return
        self.molgridview.on_mol_double_click(citem)

    def clear_canvas(self):
        self.filenames = []
        self.molgridview.clear()
        self.exit_file()
        self.statusBar().showMessage("Molecules removed", 2000)

    def closeEvent(self, event):
        self.molgridview.clear()
        self.exit_file()
        event.ignore()

    def exit_file(self):

        # self.substructure_selector.close()
        exit(0)  #TODO, how to exit QtWidgets.QApplication from within class instance?

    def about_help(self):
        QtWidgets.QMessageBox.about(
            self, "About maligner",
            """maligner is an Open-Source Molecular Alignment Tool.\n\n\nBased on RDKit: http://www.rdkit.org/\nBased on rdeditor: https://github.com/EBjerrum/rdeditor\nSome icons from: http://icons8.com\nSource code: https://github.com/hellmrf/maligner\n\nReleased under GPL-v3.0."""
        )

    def set_anchor(self):
        item = self.molgridview.listview.currentItem()
        if item is None:
            self.statusBar().showMessage("No molecule selected", 2000)
            return

        crow = self.molgridview.listview.currentRow()
        self.molgridview.set_anchor(crow)
        self.molgridview.populate_listwidget()

    def openSubsSelector(self):
        pass
        # self.substructure_selector.show()

    def set_log_level(self):
        pass
        # loglevel = self.sender().objectName().split(':')[-1].upper()
        # self.editor.logger.setLevel(loglevel)

    # Function to create actions for menus and toolbars
    def create_actions(self):
        self.openAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap("open.png")),
                                        'O&pen',
                                        self,
                                        shortcut=QtGui.QKeySequence.Open,
                                        statusTip="Open an existing file",
                                        triggered=self.open_file)

        self.saveAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap("icons8-Save.png")),
                                        'S&ave',
                                        self,
                                        shortcut=QtGui.QKeySequence.Save,
                                        statusTip="Save file",
                                        triggered=self.saveFile)

        self.exitAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap("icons8-Shutdown.png")),
                                        'E&xit',
                                        self,
                                        shortcut="Ctrl+Q",
                                        statusTip="Exit the Application",
                                        triggered=self.exit_file)

        self.aboutAction = QtGui.QAction(QtGui.QIcon(self.get_pixmap("about.png")),
                                         'A&bout',
                                         self,
                                         statusTip="Displays info about text editor",
                                         triggered=self.about_help)

        self.aboutQtAction = QtGui.QAction("About &Qt",
                                           self,
                                           statusTip="Show the Qt library's About box",
                                           triggered=QtWidgets.QApplication.aboutQt)

        #Misc Actions
        self.deleteMoleculeAction = QtGui.QAction(
            QtGui.QIcon(self.get_pixmap("icons8-Trash.png")),
            'Delete &X',
            self,
            shortcut="Ctrl+X",
            statusTip="Remove this molecule from canvas Ctrl+X",
            triggered=self.clear_canvas,
            objectName="Clear Canvas")

        self.anchorAction = QtGui.QAction(
            QtGui.QIcon(self.get_pixmap('icons8-Anchor.png')),
            'Anchor current molecule &A',
            self,
            shortcut="A",
            statusTip="Set the selected molecule as the anchor for the alignment. (A)",
            triggered=self.set_anchor,
            objectName="Set Anchor")

        self.openSelectorAction = QtGui.QAction(
            QtGui.QIcon(self.get_pixmap('icons8-Molecule.png')),
            'Open Selector',
            self,
            shortcut="S",
            statusTip="Opens the molecule selector for some molecule",
            triggered=self.open_selector,
            objectName="Open selector")

        self.loglevelactions = {}
        for key in self.loglevels:
            self.loglevelactions[key] = QtGui.QAction(key,
                                                      self,
                                                      statusTip="Set logging level to %s" % key,
                                                      triggered=self.set_log_level,
                                                      objectName="loglevel:%s" % key)
        self.loglevelactions["Debug"].setChecked(True)


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    myApp = QtWidgets.QApplication(sys.argv)
    mainWindow = MainWindow(loglevel=loglevel)
    sys.exit(myApp.exec())


if __name__ == '__main__':
    launch(loglevel="DEBUG")
