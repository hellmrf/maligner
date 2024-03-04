#!/usr/bin/env python
from __future__ import print_function
from pathlib import Path

# Import required modules
import sys, os
from typing import Optional
from PySide6.QtGui import *
from PySide6.QtWidgets import *
from PySide6 import QtWidgets

#Import model
from maligner.molEditWidget import MolEditWidget
from maligner.substructure_selector import SubstructureSelectorDialog

from rdkit import Chem


# The main window class
class MainWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self,
                 filename: Optional[Path | str] = None,
                 loglevel="WARNING"):
        super(MainWindow, self).__init__()
        self.pixmappath = os.path.abspath(
            os.path.dirname(__file__)) + '/pixmaps/'
        self.loglevels = [
            "Critical", "Error", "Warning", "Info", "Debug", "Notset"
        ]
        self.editor = MolEditWidget()
        self.substructure_selector = SubstructureSelectorDialog(
            filename=filename)
        self._filename = filename
        self.initGUI()
        # TODO: selectionChanges ainda não existe
        # self.substructure_selector.selectionChanged.connect(self.setAtomTypeName)
        self.editor.logger.setLevel(loglevel)

    #Properties
    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, filename):
        if filename != self._filename:
            self._filename = filename
            self.setWindowTitle(str(filename))

    def initGUI(self):
        self.setWindowTitle("maligner - An Open-Source Molecular Alignment Tool")
        self.setWindowIcon(QIcon(self.pixmappath + 'appicon.svg.png'))
        self.setGeometry(100, 100, 200, 150)

        # self.center = self.editor
        # self.center.setFixedSize(600, 600)
        # self.setCentralWidget(self.center)

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.SetupComponents()

        self.infobar = QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        if self.filename is not None:
            self.editor.logger.info("Loading model from %s" % self.filename)
            self.loadMolFile(self.filename)

        self.editor.sanitizeSignal.connect(self.infobar.setText)
        self.show()

    # Function to setup status bar, central widget, menu bar, tool bar
    def SetupComponents(self):
        self.myStatusBar = QStatusBar()
        #        self.molcounter = QLabel("-/-")
        #        self.myStatusBar.addPermanentWidget(self.molcounter, 0)
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)

        self.CreateActions()
        self.CreateMenus()
        self.CreateToolBars()

    # Actual menu bar item creation
    def CreateMenus(self):
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

    def CreateToolBars(self):
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

    def loadMolFile(self, filename):
        self.filename = filename
        mol = Chem.MolFromMolFile(str(self.filename),
                                  sanitize=False,
                                  strictParsing=False)
        self.editor.mol = mol
        self.statusBar().showMessage("File opened")

    def openFile(self):
        self.filename, self.filterName = QFileDialog.getOpenFileName(
            self, caption="Open MOL file", filter=self.filters)
        self.loadMolFile(self.filename)

    def saveFile(self):
        if self.filename != None:
            Chem.MolToMolFile(self.editor.mol, str(self.filename))
        else:
            self.saveAsFile()

    def saveAsFile(self):
        self.filename, self.filterName = QFileDialog.getSaveFileName(
            self, filter=self.filters)
        if (self.filename != ''):
            if self.filename[-4:].upper() != ".MOL":
                self.filename = self.filename + ".mol"
            Chem.MolToMolFile(self.editor.mol, str(self.filename))
            #            file = open(self.fileName, 'w')
            #            file.write(self.textEdit.toPlainText())
            self.statusBar().showMessage("File saved", 2000)

    def open_selector(self):
        self.substructure_selector.show()

    def clearCanvas(self):
        self.editor.mol = None
        self.filename = None
        self.exitFile()
        self.statusBar().showMessage("Molecule removed from canvas", 2000)

    def closeEvent(self, event):
        self.exitFile()
        event.ignore()

    def exitFile(self):
        self.substructure_selector.close()
        exit(0)  #TODO, how to exit qapplication from within class instance?

    def aboutHelp(self):
        QMessageBox.about(
            self, "About maligner",
            """maligner is an Open-Source Molecular Alignment Tool.\n\n\nBased on RDKit: http://www.rdkit.org/\nBased on rdeditor: https://github.com/EBjerrum/rdeditor\nSome icons from: http://icons8.com\nSource code: https://github.com/hellmrf/maligner\n\nReleased under GPL-v3.0."""
        )

    def setAction(self):
        sender = self.sender()
        self.editor.setAction(sender.objectName())
        self.myStatusBar.showMessage("Action %s selected" %
                                     sender.objectName())

    def setBondType(self):
        sender = self.sender()
        self.editor.setBondType(sender.objectName())
        self.myStatusBar.showMessage("Bondtype %s selected" %
                                     sender.objectName())

    def setAtomType(self):
        sender = self.sender()
        self.editor.setAtomType(sender.objectName())
        self.myStatusBar.showMessage("Atomtype %s selected" %
                                     sender.objectName())

    def setAtomTypeName(self, atomname):
        self.editor.setAtomType(str(atomname))
        self.myStatusBar.showMessage("Atomtype %s selected" % atomname)

    def set_anchor(self):
        raise NotImplementedError("")

    def openSubsSelector(self):
        self.substructure_selector.show()

    def setLogLevel(self):
        loglevel = self.sender().objectName().split(':')[-1].upper()
        self.editor.logger.setLevel(loglevel)

    # Function to create actions for menus and toolbars
    def CreateActions(self):
        self.openAction = QAction(QIcon(self.pixmappath + 'open.png'),
                                  'O&pen',
                                  self,
                                  shortcut=QKeySequence.Open,
                                  statusTip="Open an existing file",
                                  triggered=self.openFile)

        self.saveAction = QAction(QIcon(self.pixmappath + '/icons8-Save.png'),
                                  'S&ave',
                                  self,
                                  shortcut=QKeySequence.Save,
                                  statusTip="Save file",
                                  triggered=self.saveFile)

        self.exitAction = QAction(QIcon(self.pixmappath +
                                        'icons8-Shutdown.png'),
                                  'E&xit',
                                  self,
                                  shortcut="Ctrl+Q",
                                  statusTip="Exit the Application",
                                  triggered=self.exitFile)

        self.aboutAction = QAction(QIcon(self.pixmappath + 'about.png'),
                                   'A&bout',
                                   self,
                                   statusTip="Displays info about text editor",
                                   triggered=self.aboutHelp)

        self.aboutQtAction = QAction(
            "About &Qt",
            self,
            statusTip="Show the Qt library's About box",
            triggered=QApplication.aboutQt)

        #Misc Actions
        self.deleteMoleculeAction = QAction(
            QIcon(self.pixmappath + 'icons8-Trash.png'),
            'Delete &X',
            self,
            shortcut="Ctrl+X",
            statusTip="Remove this molecule from canvas Ctrl+X",
            triggered=self.clearCanvas,
            objectName="Clear Canvas")

        self.anchorAction = QAction(
            QIcon(self.pixmappath + 'icons8-Anchor.png'),
            'Anchor current molecule &A',
            self,
            shortcut="A",
            statusTip="Set the selected molecule as the anchor for the alignment. (A)",
            triggered=self.set_anchor,
            objectName="Set Anchor")

        self.openSelectorAction = QAction(
            QIcon(self.pixmappath + 'icons8-Molecule.png'),
            'Open Selector',
            self,
            shortcut="S",
            statusTip="Opens the molecule selector for some molecule",
            triggered=self.open_selector,
            objectName="Open selector")

        self.loglevelactions = {}
        for key in self.loglevels:
            self.loglevelactions[key] = QAction(
                key,
                self,
                statusTip="Set logging level to %s" % key,
                triggered=self.setLogLevel,
                objectName="loglevel:%s" % key)
        self.loglevelactions["Debug"].setChecked(True)


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    myApp = QApplication(sys.argv)
    argv1 = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if argv1 and argv1.is_file():
        mainWindow = MainWindow(filename=argv1, loglevel=loglevel)
    else:
        mainWindow = MainWindow(filename=None, loglevel=loglevel)
    sys.exit(myApp.exec())


if __name__ == '__main__':
    launch(loglevel="DEBUG")
