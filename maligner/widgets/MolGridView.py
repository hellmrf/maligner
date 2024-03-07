# Just copied from https://stackoverflow.com/a/62031429/5160230, which is an image gallery.
# We can just adapt this for our needs.

import os
from PySide6 import QtCore, QtGui, QtWidgets

ICON_SIZE = 100


class StyledItemDelegate(QtWidgets.QStyledItemDelegate):

    def initStyleOption(self, option, index):
        super().initStyleOption(option, index)
        option.text = option.fontMetrics.elidedText(index.data(), QtCore.Qt.ElideRight, ICON_SIZE)


class MolGridViewWidget(QtWidgets.QMainWindow):

    def __init__(self, parent=None):
        super().__init__(parent)

        self.choose_btn = QtWidgets.QPushButton(self.tr("Choose"),
                                                clicked=self.on_choose_btn_clicked)
        self.choose_btn.setFixedSize(50, 50)
        self.path_le = QtWidgets.QLineEdit()
        self.back_btn = QtWidgets.QPushButton(self.tr("^"), clicked=self.on_back_btn_clicked)
        self.back_btn.setFixedSize(20, 20)
        self.pixmap_lw = QtWidgets.QListWidget(
            viewMode=QtWidgets.QListView.IconMode,
            iconSize=ICON_SIZE * QtCore.QSize(1, 1),
            movement=QtWidgets.QListView.Static,
            resizeMode=QtWidgets.QListView.Adjust,
        )
        delegate = StyledItemDelegate(self.pixmap_lw)
        self.pixmap_lw.setItemDelegate(delegate)

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        grid_layout = QtWidgets.QGridLayout(central_widget)

        grid_layout.addWidget(self.choose_btn, 0, 0, 2, 1)
        grid_layout.addWidget(self.path_le, 0, 1)
        grid_layout.addWidget(self.back_btn, 1, 1, alignment=QtCore.Qt.AlignRight)
        grid_layout.addWidget(self.pixmap_lw, 2, 0, 1, 2)

        self.resize(640, 480)

        self.timer_loading = QtCore.QTimer(interval=50, timeout=self.load_image)
        self.filenames_iterator = None

    @QtCore.Slot()
    def on_choose_btn_clicked(self):

        directory = QtWidgets.QFileDialog.getExistingDirectory(
            options=QtWidgets.QFileDialog.DontUseNativeDialog)
        if directory:
            self.start_loading(directory)

    @QtCore.Slot()
    def on_back_btn_clicked(self):
        directory = os.path.dirname(self.path_le.text())
        self.start_loading(directory)

    def start_loading(self, directory):
        if self.timer_loading.isActive():
            self.timer_loading.stop()
        self.path_le.setText(directory)
        self.filenames_iterator = self.load_images(directory)
        self.pixmap_lw.clear()
        self.timer_loading.start()

    @QtCore.Slot()
    def load_image(self):
        try:
            filename = next(self.filenames_iterator)
        except StopIteration:
            self.timer_loading.stop()
        else:
            name = os.path.basename(filename)
            it = QtWidgets.QListWidgetItem(name)
            it.setIcon(QtGui.QIcon(filename))
            self.pixmap_lw.addItem(it)

    def load_images(self, directory):
        it = QtCore.QDirIterator(
            directory,
            ["*.jpg", "*.png"],
            QtCore.QDir.Files,
            QtCore.QDirIterator.Subdirectories,
        )
        while it.hasNext():
            filename = it.next()
            yield filename


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    w = MolGridViewWidget()
    w.show()
    sys.exit(app.exec())
