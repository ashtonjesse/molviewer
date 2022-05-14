import requests
import nglview
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Union, List, Dict
from nglview.component import ComponentViewer
from pathlib import Path

class Macromolecule():
    """
    A macromolecule obtained from the WorldWide Protein Data Bank (searchable at www.rcsb.org)
    """
    PDB_FILE_DOWNLOAD_URL = "https://files.rcsb.org/download/"

    def __init__(self, pdbid: str) -> None:
        self.pdbid = pdbid

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None) -> nglview.NGLWidget:
        if existing_viewer == None:
            return nglview.show_pdbid(self.pdbid)
        else:
            existing_viewer.add_component(nglview.adaptor.PdbIdStructure(self.pdbid))
            return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        request_response = requests.get(f"{self.PDB_FILE_DOWNLOAD_URL}{self.pdbid}.pdb")
        with open(file_path, "w") as f:
            f.write(str(request_response.content))
        return Path(file_path)


class ChemicalMolecule():
    """
    A small chemical molecule from the ChemBL database (searchable at www.ebi.ac.uk/chembl
    """
    def __init__(self, chemblid: str) -> None:
        self.chemblid = chemblid
        self.compounds_api = self.get_compounds_api()
        self.mol_data = self.get_mol_data(chemblid)
        self.conformer = self.get_conformer()

    def get_compounds_api(self) -> new_client.molecule:
        return new_client.molecule

    def get_mol_data(self, chemblid: str) -> List[Dict]:
        return self.compounds_api.filter(chembl_id=chemblid).only(['molecule_chembl_id', 'molecule_structures'])

    def get_conformer(self) -> Chem.rdchem.Mol:
        mol = Chem.AddHs(Chem.MolFromSmiles(self.mol_data[0]['molecule_structures']['canonical_smiles']))
        _ = AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        return mol

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None) -> nglview.NGLWidget:
        if existing_viewer == None:
            return nglview.show_rdkit(self.conformer)
        else:
            existing_viewer.add_component(nglview.adaptor.RdkitStructure(self.conformer))
            return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        with Chem.SDWriter(file_path) as writer:
            for cid in range(self.conformer.GetNumConformers()):
                writer.write(self.conformer, confId=cid)
        return Path(file_path)