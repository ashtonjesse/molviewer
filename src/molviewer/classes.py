import requests
import chembl_webresource_client as chembl
import rdkit_pypi as rdkit
import biopython
from typing import Union
from nglview.component import ComponentViewer
from pathlib import Path


class Molecule:
    def __init__(self) -> None:
        pass

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        pass


class Macromolecule(Molecule):
    def __init__(self, pbdid: str) -> None:
        super().__init__()
        self.pdbid = pbdid

    def save(self, file_path: Union[str, Path]) -> Path:
        pass


class ChemicalMolecule(Molecule):
    def __init__(self, chemblid: str) -> None:
        super().__init__()
        self.chemblid = chemblid

    def save(self, file_path: Union[str, Path]) -> Path:
        pass





