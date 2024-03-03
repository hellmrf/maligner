
from malign.mainwindow import MainWindow
import sys, os
from PySide6.QtWidgets import QApplication

def main(loglevel="WARNING"):
    """Function that launches the mainWindow Application"""
    try:
        app = QApplication(sys.argv)
        main_window = MainWindow(loglevel="DEBUG")
        main_window.show()
        sys.exit(app.exec())
    except Exception as e:
        print(sys.exc_info()[1])
        raise e

if __name__ == '__main__':
    main(loglevel="DEBUG")
