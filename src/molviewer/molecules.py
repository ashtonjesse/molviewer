""" Molecule classes.

Classes:
--------
Macromolecule: a molecule from the www.rcsb.org database.
Chemicalmolecule: a molecule from the www.ebi.ac.uk/chembl database.

"""
from typing import Union
from pathlib import Path

import chembl_webresource_client.query_set
from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
from chembl_webresource_client.new_client import new_client  # type: ignore
import requests
import nglview  # type: ignore


class Macromolecule():

    """ A macromolecule obtained from the WorldWide Protein Data Bank (
    searchable at www.rcsb.org).

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

    PDB_FILE_URL = "https://files.rcsb.org/download/"  # type: str

    def __init__(self, pdbid: str) -> None:
        """ Construct an object of Macromolecule class.

        :param pdbid:
            the pdbid of the corresponding macromolecule in the
            WorldWide Protein Data Bank

        """
        self.pdbid = pdbid

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None
             ) -> nglview.NGLWidget:
        """ Display the Macromolecule 3D structure inside a nglview.NGLWidget
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
        if existing_viewer is None:
            existing_viewer = nglview.show_pdbid(self.pdbid)
        elif (existing_viewer is not None and
                isinstance(existing_viewer, nglview.NGLWidget)):
            existing_viewer.add_component(
                nglview.adaptor.PdbIdStructure(self.pdbid))
        else:
            existing_viewer = None
        return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        """ Save the Macromolecule using the .pdb file format.

        :param file_path:
            The path to be used when saving the file.
        :return: Path(file_path):
            The path to the saved file

        """
        request_response = requests.get(
            f"{self.PDB_FILE_URL}{self.pdbid}.pdb")  # type: requests.Response
        try:
            with open(file_path, "w", encoding='utf-8') as file:
                file.write(str(request_response.content))
        except (FileNotFoundError, OSError):
            print("File Error: File path does not exist.")
        return Path(file_path)


class ChemicalMolecule():

    """ A small chemical molecule from the ChemBL database (searchable at
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
        """ Construct an object of Chemicalmolecule class.

        :param chemblid:
            The ChEMBLid of the corresponding small chemical molecule in
            the CHEMBL database.

        """
        self.chemblid = chemblid
        self.mol_data = self.__get_mol_data()
        self.conformer = self.__get_conformer()

    def __get_mol_data(self) -> chembl_webresource_client.query_set.QuerySet:
        """ Get the Chemicalmolecule structure from the ChEMBL database
        using the chembl_webresource_client.new_client.molecule API.

        :return:
            A List[Dict] containing the 'molecule_chembl_id' and
            'molecule_structures' information.

        """

        return new_client.molecule.filter(chembl_id=self.chemblid).only(
            ['molecule_chembl_id', 'molecule_structures'])

    def __get_conformer(self) -> Chem.rdchem.Mol:
        """ Get the 3D conformer molecule structure from the 2D
        structure specified in the 'canonical_smiles' field of self.mol_data.

        :return:
            A 3D conformer

        """
        mol = Chem.AddHs(
            Chem.MolFromSmiles(
                self.mol_data[0]['molecule_structures']['canonical_smiles'])
            )  # type: Chem.rdchem.Mol
        _ = AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True,
                                       useBasicKnowledge=True)
        return mol

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None
             ) -> nglview.NGLWidget:
        """ Display the Chemicalmolecule 3D structure inside a
        nglview.NGLWidget viewer in a Jupyter notebook.

        If no viewer currently exists, a new one is created and the 3D
        chemical molecule structure is displayed. If a viewer does currently
        exist then the 3D molecule structure is added to the existing viewer.

        :param existing_viewer:
            None if no viewer has previously been created or an
            nglview.NGLWidget viewer.
        :return: existing_viewer:
            An nglview.NGLWidget viewer.

        """
        if existing_viewer is None:
            existing_viewer = nglview.show_rdkit(
                self.conformer)
        elif (existing_viewer is not None and
              isinstance(existing_viewer, nglview.NGLWidget)):
            existing_viewer.add_component(
                nglview.adaptor.RdkitStructure(self.conformer))
        else:
            existing_viewer = None
        return existing_viewer

    def save(self, file_path: Union[str, Path]) -> Path:
        """ Save the Chemicalmolecule 3D structure using the .sdf file format.

        :param file_path:
            The path to be used when saving the file.
        :return: Path(file_path):
            The path to the saved file.

        """
        try:
            with Chem.SDWriter(file_path) as writer:
                for cid in range(self.conformer.GetNumConformers()):
                    writer.write(self.conformer, confId=cid)
        except (FileNotFoundError, OSError):
            print("File error: File path does not exist.")
        return Path(file_path)
