import requests
from chembl_webresource_client.new_client import new_client
from rdkit.Chem import AllChem as Chem
import Bio as biopython
import json
import jupyterlab

from typing import Union, List, Dict
from nglview.component import ComponentViewer
from pathlib import Path


class Molecule:
    def __init__(self) -> None:
        pass




class Macromolecule(Molecule):
    """
    A macromolecule obtained from the WorldWide Protein Data Bank (searchable at www.rcsb.org)
    """
    def __init__(self, pdbid: str) -> None:
        super().__init__()
        self.pdbid = pdbid

    def save(self, file_path: Union[str, Path]) -> Path:
        pass

class ChemicalMolecule(Molecule):
    """
    A small chemical molecule from the ChemBL database (searchable at www.ebi.ac.uk/chembl
    """
    def __init__(self, chemblid: str) -> None:
        super().__init__()
        self.chemblid = chemblid
        self.compounds_api = self.get_compounds_api()
        self.mol_data = self.get_mol_data(chemblid)
        self.conformer = self.get_conformer()

    def get_compounds_api(self) -> new_client.molecule:
        return new_client.molecule

    def get_mol_data(self, chemblid: str) -> Dict:
        return self.compounds_api.filter(chembl_id=chemblid).only(['molecule_chembl_id', 'pref_name', 'molecule_structures'])

    def get_conformer(self) -> rdkit.Chem.rdchem.Mol:
        mol = Chem.AddHs(Chem.MolFromSmiles(self.mol_data['molecular_structures']['canonical_smiles']))
        _ = AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        return mol

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        return nv.show_rdkit(self.conformer)

    def save(self, file_path: Union[str, Path]) -> Path:
        pass