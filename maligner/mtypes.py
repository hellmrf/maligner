from dataclasses import dataclass
from typing import TypeAlias
from pathlib import Path

from PySide6 import QtGui
from rdkit import Chem

Mol: TypeAlias = Chem.rdchem.Mol


@dataclass
class MolData:
    mol: Mol
    name: str
    filename: Path
    qicon: QtGui.QIcon
    anchor: bool = False
