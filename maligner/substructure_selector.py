#!/usr/bin/env python
from __future__ import print_function
from pathlib import Path

# Import required modules
import sys, os
from typing import Optional
from PySide6.QtGui import *
from PySide6.QtWidgets import *
from PySide6 import QtWidgets

from PySide6.QtWidgets import QDialog, QVBoxLayout

#Import model
from maligner.molEditWidget import MolEditWidget

from rdkit import Chem


# The main window class
class SubstructureSelectorWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, filename: Optional[Path | str], loglevel="WARNING"):
        super(SubstructureSelectorWindow, self).__init__()
        self.pixmappath = os.path.abspath(
            os.path.dirname(__file__)) + '/pixmaps/'
        self.loglevels = [
            "Critical", "Error", "Warning", "Info", "Debug", "Notset"
        ]
        self.editor = MolEditWidget()
        self._filename = filename
        self.initGUI(filename=filename)
        # TODO: conectar signal selectionChanged
        # self.editor.selectionChanged.connect(self.setAtomTypeName)
        self.editor.logger.setLevel(loglevel)

    #Properties
    @property
    def fileName(self):
        return self._filename

    @fileName.setter
    def fileName(self, filename):
        if filename != self._filename:
            self._filename = filename
            self.setWindowTitle(str(filename))

    def initGUI(self, filename=None):
        self.setWindowTitle("Substructure selector")
        self.setWindowIcon(QIcon(self.pixmappath + 'icons8-Cursor.png'))
        self.setGeometry(100, 100, 200, 150)

        self.center = self.editor
        self.center.setFixedSize(600, 600)
        self.setCentralWidget(self.center)
        self.fileName = filename

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.SetupComponents()

        self.infobar = QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        if self.fileName is not None:
            self.editor.logger.info("Loading model from %s" % self.fileName)
            self.loadMolFile(filename)

        self.editor.sanitizeSignal.connect(self.infobar.setText)
        self.show()

    # Function to setup status bar, central widget, menu bar, tool bar
    def SetupComponents(self):
        self.myStatusBar = QStatusBar()
        self.setStatusBar(self.myStatusBar)
        self.myStatusBar.showMessage('Ready', 10000)

        self.CreateActions()
        self.CreateToolBars()

    def CreateToolBars(self):
        self.mainToolBar = self.addToolBar('Main')
        #Main action bar
        self.mainToolBar.addAction(self.cleanCoordinatesAction)
        self.mainToolBar.addAction(self.undoAction)
        self.mainToolBar.addAction(self.deleteMoleculeAction)

    def loadMolFile(self, filename):
        self.fileName = filename
        mol = Chem.MolFromMolFile(str(self.fileName),
                                  sanitize=False,
                                  strictParsing=False)
        self.editor.mol = mol
        self.statusBar().showMessage("File opened")

    def openFile(self):
        self.fileName, self.filterName = QFileDialog.getOpenFileName(
            self, caption="Open MOL file", filter=self.filters)
        self.loadMolFile(self.fileName)

    def saveFile(self):
        if self.fileName != None:
            Chem.MolToMolFile(self.editor.mol, str(self.fileName))
        else:
            self.saveAsFile()

    def saveAsFile(self):
        self.fileName, self.filterName = QFileDialog.getSaveFileName(
            self, filter=self.filters)
        if (self.fileName != ''):
            if self.fileName[-4:].upper() != ".MOL":
                self.fileName = self.fileName + ".mol"
            Chem.MolToMolFile(self.editor.mol, str(self.fileName))
            # file = open(self.fileName, 'w')
            # file.write(self.textEdit.toPlainText())
            self.statusBar().showMessage("File saved", 2000)

    def clearCanvas(self):
        self.editor.mol = None
        self.fileName = None
        self.exitFile()
        self.statusBar().showMessage("Molecule removed from canvas", 2000)

    def closeEvent(self, event):
        self.exitFile()
        event.ignore()

    def exitFile(self):
        exit(0)  #TODO, how to exit qapplication from within class instance?

    # Function to show Diaglog box with provided Title and Message
    def msgApp(self, title, msg):
        userInfo = QMessageBox.question(self, title, msg,
                                        QMessageBox.Yes | QMessageBox.No)
        if userInfo == QMessageBox.Yes:
            return "Y"
        if userInfo == QMessageBox.No:
            return "N"
        self.close()

    def aboutHelp(self):
        QMessageBox.about(
            self, "About Simple Molecule Editor",
            """A Simple Molecule Editor where you can edit molecules\nBased on RDKit! http://www.rdkit.org/ \nSome icons from http://icons8.com\n\nSource code: https://github.com/EBjerrum/rdeditor"""
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

    def setLogLevel(self):
        loglevel = self.sender().objectName().split(':')[-1].upper()
        self.editor.logger.setLevel(loglevel)

    # Function to create actions for menus and toolbars
    def CreateActions(self):
        self.exitAction = QAction(QIcon(self.pixmappath +
                                        'icons8-Shutdown.png'),
                                  'E&xit',
                                  self,
                                  shortcut="Ctrl+Q",
                                  statusTip="Exit the Application",
                                  triggered=self.exitFile)

        #Misc Actions
        self.undoAction = QAction(
            QIcon(self.pixmappath + 'prev.png'),
            'U&ndo',
            self,
            shortcut="Ctrl+Z",
            statusTip="Undo/Redo changes to molecule Ctrl+Z",
            triggered=self.editor.undo,
            objectName="undo")

        self.deleteMoleculeAction = QAction(
            QIcon(self.pixmappath + 'icons8-Trash.png'),
            'Delete molecule from set (Del)',
            self,
            shortcut="Del",
            statusTip="Remove this molecule from set (Del)",
            triggered=self.clearCanvas,
            objectName="Clear Canvas")

        self.cleanCoordinatesAction = QAction(
            QIcon(self.pixmappath + 'icons8-Broom.png'),
            'Recalculate coordinates (F9)',
            self,
            shortcut="F9",
            statusTip="Re-calculates coordinates and redraw (F9)",
            triggered=self.editor.canon_coords_and_draw,
            objectName="Recalculate Coordinates")


class SubstructureSelectorDialog(QDialog):

    def __init__(self, filename: Optional[Path | str], loglevel="WARNING"):
        super().__init__()
        self.ss_window = SubstructureSelectorWindow(
            filename, loglevel)
        layout = QVBoxLayout()
        layout.addWidget(self.ss_window)
        self.setLayout(layout)
        self.setWindowTitle("Substructure Selector")
        self.setWindowIcon(QIcon(self.ss_window.pixmappath + 'icons8-Cursor.png'))


def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    myApp = QApplication(sys.argv)
    argv1 = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if argv1 and argv1.is_file():
        mainWindow = SubstructureSelectorWindow(filename=argv1,
                                                loglevel=loglevel)
    else:
        mainWindow = SubstructureSelectorWindow(filename=None,
                                                loglevel=loglevel)
    sys.exit(myApp.exec())


if __name__ == '__main__':
    launch(loglevel="DEBUG")
