from dataclasses import dataclass, field
from pathlib import Path
from typing import TypeAlias

from PySide6 import QtGui
from rdkit import Chem

Mol: TypeAlias = Chem.Mol


@dataclass
class MolData:
    mol: Mol
    name: str
    filename: Path
    qicon: QtGui.QIcon
    anchor: bool = False
    selected: list[int] = field(default_factory=list)
