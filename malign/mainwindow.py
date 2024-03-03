#!/usr/bin/env python
from __future__ import print_function


# Import required modules
import sys, time, os
from PySide6.QtGui import *
from PySide6.QtWidgets import *
from PySide6.QtCore import QByteArray
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6 import QtSvg

#Import model
from malign.molEditWidget import MolEditWidget
from malign.ptable_widget import PTable

from rdkit import Chem

# The main window class
class MainWindow(QtWidgets.QMainWindow):
    # Constructor function
    def __init__(self, fileName=None, loglevel="WARNING"):
        super(MainWindow,self).__init__()
        self.pixmappath = os.path.abspath(os.path.dirname(__file__)) + '/pixmaps/'
        self.loglevels = ["Critical","Error","Warning","Info","Debug","Notset"]
        self.editor = MolEditWidget()
        self.ptable = PTable()
        self._fileName = None
        self.initGUI(fileName = fileName)
        self.ptable.atomtypeChanged.connect(self.setAtomTypeName)
        self.editor.logger.setLevel(loglevel)

    #Properties
    @property
    def fileName(self):
        return self._fileName

    @fileName.setter
    def fileName(self, filename):
        if filename != self._fileName:
            self._fileName = filename
            self.setWindowTitle(str(filename))


    def initGUI(self, fileName=None):
        self.setWindowTitle("A simple mol editor")
        self.setWindowIcon(QIcon(self.pixmappath + 'appicon.svg.png'))
        self.setGeometry(100, 100, 200, 150)

        self.center = self.editor
        self.center.setFixedSize(600,600)
        self.setCentralWidget(self.center)
        self.fileName = fileName

        self.filters = "MOL Files (*.mol *.mol);;Any File (*)"
        self.SetupComponents()

        self.infobar = QLabel("")
        self.myStatusBar.addPermanentWidget(self.infobar, 0)

        if self.fileName is not None:
            self.editor.logger.info("Loading model from %s"%self.fileName)
            self.loadMolFile(fileName)

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
        self.atomtypeMenu = self.menuBar().addMenu("&AtomTypes")
        self.bondtypeMenu = self.menuBar().addMenu("&BondTypes")
        self.helpMenu = self.menuBar().addMenu("&Help")

        # File
        self.fileMenu.addAction(self.openAction)
        self.fileMenu.addAction(self.saveAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        # Tools
        self.toolMenu.addSeparator()
        self.toolMenu.addAction(self.cleanCoordinatesAction)
        self.toolMenu.addSeparator()
        self.toolMenu.addAction(self.undoAction)
        self.toolMenu.addSeparator()
        # self.toolMenu.addAction(self.removeAction)

        #Help menu
        self.helpMenu.addAction(self.aboutAction)
        self.helpMenu.addSeparator()
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
        self.mainToolBar.addAction(self.cleanCoordinatesAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.deleteMoleculeAction)
        self.mainToolBar.addSeparator()
        self.mainToolBar.addAction(self.undoAction)

    def loadMolFile(self, filename):
        self.fileName = filename
        mol = Chem.MolFromMolFile(str(self.fileName), sanitize=False, strictParsing=False)
        self.editor.mol = mol
        self.statusBar().showMessage("File opened")

    def openFile(self):
        self.fileName, self.filterName = QFileDialog.getOpenFileName(self, caption = "Open MOL file",filter = self.filters)
        self.loadMolFile(self.fileName)

    def saveFile(self):
        if self.fileName != None:
            Chem.MolToMolFile(self.editor.mol, str(self.fileName))
        else:
            self.saveAsFile()

    def saveAsFile(self):
        self.fileName, self.filterName = QFileDialog.getSaveFileName(self, filter=self.filters)
        if(self.fileName != ''):
            if self.fileName[-4:].upper() != ".MOL":
                self.fileName = self.fileName + ".mol"
            Chem.MolToMolFile(self.editor.mol, str(self.fileName))
#            file = open(self.fileName, 'w')
#            file.write(self.textEdit.toPlainText())
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
        self.ptable.close()
        exit(0) #TODO, how to exit qapplication from within class instance?


    # Function to show Diaglog box with provided Title and Message
    def msgApp(self,title,msg):
        userInfo = QMessageBox.question(self,title,msg,
                                        QMessageBox.Yes | QMessageBox.No)
        if userInfo == QMessageBox.Yes:
            return "Y"
        if userInfo == QMessageBox.No:
            return "N"
        self.close()

    def aboutHelp(self):
        QMessageBox.about(self, "About Simple Molecule Editor",
                """A Simple Molecule Editor where you can edit molecules\nBased on RDKit! http://www.rdkit.org/ \nSome icons from http://icons8.com\n\nSource code: https://github.com/EBjerrum/rdeditor""")

    def setAction(self):
        sender = self.sender()
        self.editor.setAction(sender.objectName())
        self.myStatusBar.showMessage("Action %s selected"%sender.objectName())

    def setBondType(self):
        sender = self.sender()
        self.editor.setBondType(sender.objectName())
        self.myStatusBar.showMessage("Bondtype %s selected"%sender.objectName())

    def setAtomType(self):
        sender = self.sender()
        self.editor.setAtomType(sender.objectName())
        self.myStatusBar.showMessage("Atomtype %s selected"%sender.objectName())

    def setAtomTypeName(self, atomname):
        self.editor.setAtomType(str(atomname))
        self.myStatusBar.showMessage("Atomtype %s selected"%atomname)

    def openPtable(self):
        self.ptable.show()

    def setLogLevel(self):
        loglevel = self.sender().objectName().split(':')[-1].upper()
        self.editor.logger.setLevel(loglevel)



    # Function to create actions for menus and toolbars
    def CreateActions(self):
        self.openAction = QAction( QIcon(self.pixmappath + 'open.png'), 'O&pen',
                                  self, shortcut=QKeySequence.Open,
                                  statusTip="Open an existing file",
                                  triggered=self.openFile)

        self.saveAction = QAction( QIcon(self.pixmappath + '/icons8-Save.png'), 'S&ave',
                                  self, 
                                  shortcut=QKeySequence.Save,
                                  statusTip="Save file",
                                  triggered=self.saveFile)

        self.exitAction = QAction( QIcon(self.pixmappath + 'icons8-Shutdown.png'), 'E&xit',
                                   self, shortcut="Ctrl+Q",
                                   statusTip="Exit the Application",
                                   triggered=self.exitFile)

        self.aboutAction = QAction( QIcon(self.pixmappath + 'about.png'), 'A&bout',
                                    self, statusTip="Displays info about text editor",
                                   triggered=self.aboutHelp)

        self.aboutQtAction = QAction("About &Qt", self,
                                statusTip="Show the Qt library's About box",
                                triggered=QApplication.aboutQt)


        #Misc Actions
        self.undoAction = QAction( QIcon(self.pixmappath + 'prev.png'), 'U&ndo',
                           self, shortcut="Ctrl+Z",
                           statusTip="Undo/Redo changes to molecule Ctrl+Z",
                           triggered=self.editor.undo, objectName="undo")

        self.deleteMoleculeAction = QAction( QIcon(self.pixmappath + 'icons8-Trash.png'), 'Delete &X',
                                   self, shortcut="Ctrl+X",
                                   statusTip="Remove this molecule from canvas Ctrl+X",
                                   triggered=self.clearCanvas, objectName="Clear Canvas")

        self.cleanCoordinatesAction = QAction( QIcon(self.pixmappath + 'icons8-Broom.png'), 'Recalculate coordinates &F',
                                   self, shortcut="Ctrl+F",
                                   statusTip="Re-calculates coordinates and redraw",
                                   triggered=self.editor.canon_coords_and_draw, objectName="Recalculate Coordinates")

        self.loglevelactions = {}
        for key in self.loglevels:
            self.loglevelactions[key] = QAction(key,
                                   self,
                                   statusTip="Set logging level to %s"%key,
                                   triggered=self.setLogLevel, objectName="loglevel:%s"%key)
        self.loglevelactions["Debug"].setChecked(True)

def launch(loglevel="WARNING"):
    "Function that launches the mainWindow Application"
    # Exception Handling
    try:
        myApp = QApplication(sys.argv)
        try:
            mainWindow = MainWindow(fileName = sys.argv[1], loglevel=loglevel)
        except:
            mainWindow = MainWindow(loglevel=loglevel)
        myApp.exec()
        sys.exit(0)
    except NameError:
        print("Name Error:", sys.exc_info()[1])
    except SystemExit:
        print("Closing Window...")
    except Exception:
        print(sys.exc_info()[1])


if __name__ == '__main__':
    launch(loglevel="DEBUG")

