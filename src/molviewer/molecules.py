from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import AllChem
from Bio.PDB import PDBList

from typing import Union, List, Dict
from nglview.component import ComponentViewer
import nglview as nv
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

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        if existing_viewer == None:
            #existing_viewer = nv.NGLWidget()
            return nv.show_pdbid(self.pdbid)
        else:
            return existing_viewer.add_component(nv.adaptor.PdbIdStructure(self.pdbid))

    def save(self, file_path: Union[str, Path]) -> Path:
        pdbl = PDBList()
        file_out = pdbl.retrieve_pdb_file(pdb_code=self.pdbid, pdir=file_path, file_format='pdb')
        return Path(file_out)



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

    def get_mol_data(self, chemblid: str) -> List[Dict]:
        return self.compounds_api.filter(chembl_id=chemblid).only(['molecule_chembl_id', 'molecule_structures'])

    def get_conformer(self) -> Chem.rdchem.Mol:
        mol = Chem.AddHs(Chem.MolFromSmiles(self.mol_data[0]['molecule_structures']['canonical_smiles']))
        _ = AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        return mol

    def show(self, existing_viewer: Union[ComponentViewer, None] = None) -> ComponentViewer:
        if existing_viewer == None:
            #existing_viewer = nv.NGLWidget()
            return nv.show_rdkit(self.conformer)
        else:
            return existing_viewer.add_component(nv.adaptor.RdkitStructure(self.conformer))

    def save(self, file_path: Union[str, Path]) -> Path:
        with Chem.SDWriter(file_path) as writer:
            for cid in range(self.conformer.GetNumConformers()):
                writer.write(self.conformer, confId=cid)

        return Path(file_path)