import requests
import nglview # type: ignore
from chembl_webresource_client.new_client import new_client # type: ignore
from rdkit import Chem # type: ignore
from rdkit.Chem import AllChem # type: ignore
from typing import Union, List, Dict
from pathlib import Path

class Macromolecule():

    """
    A macromolecule obtained from the WorldWide Protein Data Bank (searchable
    at www.rcsb.org).

    Attributes
    ----------
    pdbid : str
        the pdbid of the corresponding macromolecule in the WorldWide Protein
        Data Bank.

    Public Methods
    -------
    show(existing_viewer): Display the Macromolecule 3D structure.
    save(file_path): Save the Macromolecule using the .pdb file format.
    
    """
    
    PDB_FILE_DOWNLOAD_URL = "https://files.rcsb.org/download/" # type: str

    def __init__(self, pdbid: str) -> None:
        """
        Construct an object of Macromolecule class.

        :param pdbid:
            the pdbid of the corresponding macromolecule in the
            WorldWide Protein Data Bank

        """
        self.pdbid = pdbid

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None) -> nglview.NGLWidget:
        """
        Display the Macromolecule 3D structure inside a nglview.NGLWidget
        viewer in a Jupyter notebook.

        If no viewer currently exists, a new one is created and the 3D
        macromolecule structure is displayed. If a viewer does currently
        exist then the 3D molecule structure is added to the existing viewer.

        :param existing_viewer:
            None if no viewer has previously been created or an
            nglview.NGLWidget viewer.
        :return: existing_viewer:
            An nglview.NGLWidget viewer.

        """
        if existing_viewer == None:
            return nglview.show_pdbid(self.pdbid)
        elif existing_viewer is not None and isinstance(existing_viewer, nglview.NGLWidget):
            existing_viewer.add_component(nglview.adaptor.PdbIdStructure(self.pdbid))
            return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        """
        Save the Macromolecule using the .pdb file format.

        :param file_path:
            The path to be used when saving the file.
        :return: Path(file_path):
            The path to the saved file

        """
        request_response = requests.get(f"{self.PDB_FILE_DOWNLOAD_URL}{self.pdbid}.pdb") # type: requests.Response
        with open(file_path, "w") as f:
            f.write(str(request_response.content))
        return Path(file_path)


class ChemicalMolecule():

    """
    A small chemical molecule from the ChemBL database (searchable at
    www.ebi.ac.uk/chembl)

    Attributes
    ----------
    chemblid : str
        The chemblid of the corresponding small chemical molecule in the
        ChEMBL database.


    Public Methods
    -------
    show(existing_viewer): Display the Chemicalmolecule 3D structure.
    save(file_path): Save the Chemicalmolecule using the .sdf file format.

    """

    def __init__(self, chemblid: str) -> None:
        """
        Construct an object of Chemicalmolecule class.

        :param chemblid:
            The ChEMBLid of the corresponding small chemical molecule in
            the CHEMBL database.

        """
        self.chemblid = chemblid
        self.compounds_api = self.__get_compounds_api()  # type: new_client.molecule
        self.mol_data = self.__get_mol_data(chemblid)  # type: List[Dict]
        self.conformer = self.__get_conformer() # type: Chem.rdchem.Mol

    def __get_compounds_api(self) -> new_client.molecule:
        """
        Get the chembl_webresource_client.new_client.molecule API from which
        the Chemicalmolecule structure can be obtained.

        :return:
            An instance of chembl_webresource_client.new_client.molecule

        """
        return new_client.molecule

    def __get_mol_data(self, chemblid: str) -> List[Dict]:
        """
        Get the Chemicalmolecule structure from the ChEMBL database using the
        chembl_webresource_client.new_client.molecule API.

        :return:
            A List[Dict] containing the 'molecule_chembl_id' and
            'molecule_structures' information.

        """
        return self.compounds_api.filter(chembl_id=chemblid).only(['molecule_chembl_id', 'molecule_structures'])

    def __get_conformer(self) -> Chem.rdchem.Mol:
        """
        Get the 3D conformer molecule structure from the 2D
        structure specified in the 'canonical_smiles' field of self.mol_data.
                
        :return:
            A 3D conformer

        """
        mol = Chem.AddHs(Chem.MolFromSmiles(self.mol_data[0]['molecule_structures']['canonical_smiles'])) # type: Chem.rdchem.Mol
        _ = AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        return mol

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None) -> nglview.NGLWidget:
        """
        Display the Chemicalmolecule 3D structure inside a nglview.NGLWidget
        viewer in a Jupyter notebook.

        If no viewer currently exists, a new one is created and the 3D
        chemical molecule structure is displayed. If a viewer does currently
        exist then the 3D molecule structure is added to the existing viewer.

        :param existing_viewer:
            None if no viewer has previously been created or an
            nglview.NGLWidget viewer.
        :return: existing_viewer:
            An nglview.NGLWidget viewer.

        """
        if existing_viewer == None:
            return nglview.show_rdkit(self.conformer)
        elif existing_viewer is not None and isinstance(existing_viewer, nglview.NGLWidget):
            existing_viewer.add_component(nglview.adaptor.RdkitStructure(self.conformer))
            return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        """
        Save the Chemicalmolecule 3D structure using the .sdf file format.

        :param file_path:
            The path to be used when saving the file.
        :return: Path(file_path):
            The path to the saved file.

        """
        with Chem.SDWriter(file_path) as writer:
            for cid in range(self.conformer.GetNumConformers()):
                writer.write(self.conformer, confId=cid)
        return Path(file_path)