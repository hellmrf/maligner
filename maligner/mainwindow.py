import sys
from pathlib import Path
from typing import Optional

from PySide6 import QtCore, QtGui, QtWidgets

from maligner.icons import pixmap
from maligner.widgets.MolGridView import MolGridViewWidget


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, filenames: Optional[list[Path | str]] = None):
        super(MainWindow, self).__init__()
        self._filenames = filenames
        self.molgridview = MolGridViewWidget(self)
        self.init_GUI()

    def tr(self, sourceText: str, disambiguation: Optional[str] = None, n: int = -1) -> str:
        return QtWidgets.QApplication.translate("MolGridViewWidget", sourceText, disambiguation, n)

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        if filename != self._filename:
            self._filename = filename
            self.setWindowTitle(str(filename))

    def init_GUI(self):
        self.setWindowTitle(
            self.tr(r"maligner |\ An Open-Source Molecular Alignment Tool")
        )
        self.setWindowIcon(QtGui.QIcon(pixmap("appicon.svg.png")))
        self.setGeometry(400, 400, 700, 500)

        self.setCentralWidget(self.molgridview)

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.setup_components()

        self.infobar = QtWidgets.QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        self.show()

    def setup_components(self):
        self.myStatusBar = QtWidgets.QStatusBar()
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage("Ready", 10000)

        self.create_actions()
        self.create_menus()
        self.create_tool_bars()

    def create_menus(self):
        # File
        self.fileMenu = self.menuBar().addMenu(self.tr("File"))
        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addAction(self.saveAction)
        self.fileMenu.addAction(self.save_alignment_action)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        # Tools
        self.toolMenu = self.menuBar().addMenu(self.tr("Tools"))
        self.toolMenu.addAction(self.anchorAction)
        self.toolMenu.addAction(self.deleteMoleculeAction)
        self.toolMenu.addAction(self.computeMCSAction)
        self.toolMenu.addAction(self.runAlignmentAction)

        # Help menu
        self.helpMenu = self.menuBar().addMenu(self.tr("Help"))
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addAction(self.aboutQtAction)

    def create_tool_bars(self):
        self.mainToolBar = self.addToolBar("Main")
        self.mainToolBar.addAction(self.openAction)
        self.mainToolBar.addAction(self.saveAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.anchorAction)
        self.mainToolBar.addAction(self.computeMCSAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.deleteMoleculeAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.runAlignmentAction)

    def open_file(self):
        self.molgridview.file_chooser()
        self.molgridview.load_molecules()

    def open_selector(self):
        citem = self.molgridview.listview.currentItem()
        if citem is None:
            self.statusBar().showMessage(self.tr("No molecule selected"), 2000)
            return
        self.molgridview.on_mol_double_click(citem)

    def clear_canvas(self):
        self.filenames = []
        self.molgridview.clear()
        self.exit_file()
        self.statusBar().showMessage(self.tr("Molecules removed"), 2000)

    def exit_file(self):
        exit(0)

    def about_help(self):
        QtWidgets.QMessageBox.about(
            self,
            self.tr("About maligner"),
            self.tr(
                """maligner is an Open-Source Molecular Alignment Tool.\n\n\nBased on RDKit: http://www.rdkit.org/\nBased on rdeditor: https://github.com/EBjerrum/rdeditor\nSome icons from: http://icons8.com\nSource code: https://github.com/hellmrf/maligner\n\nReleased under GPL-v3.0."""
            ),
        )

    def set_anchor(self):
        item = self.molgridview.listview.currentItem()
        if item is None:
            self.statusBar().showMessage(self.tr("No molecule selected"), 2000)
            return

        crow = self.molgridview.listview.currentRow()
        self.molgridview.set_anchor(crow)
        self.molgridview.populate_listwidget()

    def compute_MCS(self):
        self.molgridview.compute_MCS()
        self.molgridview.populate_listwidget()

    def run_alignment(self):
        self.molgridview.run_alignment()
        self.molgridview.populate_listwidget()

    def save_alignment_file(self):
        self.molgridview.save_alignment_file()

    def create_actions(self):
        self.openAction = QtGui.QAction(
            QtGui.QIcon(pixmap("open.png")), self.tr("Open"), self
        )
        self.openAction.setShortcut(QtGui.QKeySequence.StandardKey.Open)
        self.openAction.setStatusTip(self.tr("Open an existing file"))
        self.openAction.triggered.connect(self.open_file)

        self.saveAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Save.png")), self.tr("Save"), self
        )
        self.saveAction.setShortcut(QtGui.QKeySequence.StandardKey.Save)
        self.saveAction.setStatusTip(self.tr("Save file"))

        self.save_alignment_action = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Save as.png")),
            self.tr("Save alignment file"),
            self,
        )
        self.save_alignment_action.setShortcut(QtGui.QKeySequence.SaveAs)
        self.save_alignment_action.setStatusTip(
            self.tr("Save the alignment to a CSV file.")
        )
        self.save_alignment_action.triggered.connect(self.save_alignment_file)

        self.exitAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Shutdown.png")), self.tr("Exit"), self
        )
        self.exitAction.setStatusTip(self.tr("Exit the Application"))
        self.exitAction.triggered.connect(self.exit_file)

        self.aboutAction = QtGui.QAction(
            QtGui.QIcon(pixmap("about.png")), self.tr("About"), self
        )
        self.aboutAction.setStatusTip(self.tr("Displays info about maligner"))
        self.aboutAction.triggered.connect(self.about_help)

        self.aboutQtAction = QtGui.QAction(self.tr("About Qt"), self)
        self.aboutQtAction.setStatusTip(self.tr("Show the Qt library's About box"))
        self.aboutQtAction.triggered.connect(QtWidgets.QApplication.aboutQt)

        self.deleteMoleculeAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Trash.png")), self.tr("Delete"), self
        )
        self.deleteMoleculeAction.setShortcut(QtGui.QKeySequence.StandardKey.Delete)
        self.deleteMoleculeAction.setStatusTip(
            self.tr("Remove this molecule from canvas")
        )
        self.deleteMoleculeAction.triggered.connect(self.clear_canvas)

        self.anchorAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Anchor.png")),
            self.tr("Anchor current molecule"),
            self,
        )
        self.anchorAction.setShortcut("A")
        self.anchorAction.setStatusTip(
            self.tr("Set the selected molecule as the anchor for the alignment. (A)")
        )
        self.anchorAction.triggered.connect(self.set_anchor)

        self.openSelectorAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Molecule.png")), "Open Selector", self
        )
        self.openSelectorAction.triggered.connect(self.open_selector)
        self.openSelectorAction.setStatusTip(
            "Opens the molecule selector for some molecule"
        )

        self.computeMCSAction = QtGui.QAction(
            QtGui.QIcon(pixmap("MCS.png")), self.tr("Compute MCS"), self
        )
        self.computeMCSAction.setShortcut("S")
        self.computeMCSAction.setStatusTip(
            self.tr(
                "Computes and select the Maximum Common Substructure for the loaded molecules."
            )
        )
        self.computeMCSAction.triggered.connect(self.compute_MCS)

        self.runAlignmentAction = QtGui.QAction(
            QtGui.QIcon(pixmap("icons8-Molecule.png")), self.tr("Run Alignment"), self
        )
        self.runAlignmentAction.setShortcut("R")
        self.runAlignmentAction.setStatusTip(self.tr("Align all the molecules."))
        self.runAlignmentAction.triggered.connect(self.run_alignment)


def launch():
    """Function that launches the mainWindow Application"""
    app = QtWidgets.QApplication(sys.argv)

    path = str(Path(__file__).parent / "translations")
    translator = QtCore.QTranslator(app)
    if translator.load(QtCore.QLocale.system(), "", "", path):
        app.installTranslator(translator)

    mainWindow = MainWindow()

    sys.exit(app.exec())


if __name__ == "__main__":
    launch()
