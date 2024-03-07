from pathlib import Path
from PySide6 import QtGui

# icons8_Cursor = Path(__file__).parent / "pixmaps" / "icons8-Cursor.png"


def pixmap(name: str) -> QtGui.QPixmap:
    p = Path(__file__).parent / 'pixmaps' / name
    return QtGui.QPixmap(p)


def icon(name: str) -> QtGui.QIcon:
    return QtGui.QIcon(pixmap(name))
