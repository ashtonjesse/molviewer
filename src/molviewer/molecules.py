"""Molviewer Molecule classes.

A molecule is either a Macromolecule with information retrieved from the
`www.rcsb.org <https://www.rcsb.org>`_ database, or a ChemicalMolecule
with information retrieved from the
`www.ebi.ac.uk/chembl <https://www.ebi.ac.uk/chembl/>`_ database.

"""

from typing import Union, ClassVar
from pathlib import Path

import chembl_webresource_client.query_set
from rdkit import Chem
from rdkit.Chem import AllChem
from chembl_webresource_client.new_client import new_client
import requests
import nglview


class Molecule:

    """The base class for molecules.

    """

    structure: Union[nglview.Structure, None] = None

    def __init__(self) -> None:
        pass

    def show(self, existing_viewer: Union[nglview.NGLWidget, None] = None
             ) -> nglview.NGLWidget:
        """Display the molecule's 3D structure inside a nglview.NGLWidget
        viewer in a Jupyter notebook. If no viewer currently exists,
        a new one is created and the 3D molecular structure is
        displayed. If a viewer does currently exist then the 3D
        molecule structure is added to the existing viewer.

        :param existing_viewer: a 3D molecule model viewer, defaults to None.
        :return: An nglview.NGLWidget viewer.

        """
        if existing_viewer is None:
            existing_viewer: nglview.NGLWidget = nglview.NGLWidget()
        if isinstance(existing_viewer,
                      nglview.NGLWidget) and self.structure is not None:
            existing_viewer.add_component(self.structure)
        return existing_viewer

    
class Macromolecule(Molecule):

    """A macromolecule obtained from the WorldWide Protein Data Bank (
    searchable at `www.rcsb.org <https://www.rcsb.org>`_).

    """

    PDB_FILE_URL: ClassVar[str] = "https://files.rcsb.org/download/"
    pdbid: str

    def __init__(self, pdbid: str) -> None:
        """Construct an object of :class:`Macromolecule`.

        :param pdbid: the pdbid of the corresponding macromolecule in the
            WorldWide Protein Data Bank

        """
        super().__init__()
        self.pdbid: str = pdbid
        self.structure: nglview.Structure = nglview.adaptor.PdbIdStructure(
            self.pdbid)

    def save(self, file_path: Union[str, Path]) -> Path:
        """Save the Macromolecule using the .pdb file format.

        :param file_path: The path to be used when saving the file.
        :return: The path to the saved file
        :raise: FileNotFoundError or OSError if the specified file path
            does not exist.

        """
        request_response: requests.Response = requests.get(
            f"{self.PDB_FILE_URL}{self.pdbid}.pdb")
        try:
            with open(file_path, "w", encoding='utf-8') as file:
                file.write(str(request_response.content))
        except (FileNotFoundError, OSError):
            print("File Error: File path does not exist.")
        return Path(file_path)


class ChemicalMolecule(Molecule):

    """A small chemical molecule from the ChemBL database (searchable at
    `www.ebi.ac.uk/chembl <https://www.ebi.ac.uk/chembl/>`_)

    """

    chemblid: str
    mol_data: chembl_webresource_client.query_set.QuerySet
    conformer: Chem.rdchem.Mol

    def __init__(self, chemblid: str) -> None:
        """ Construct an object of :class:`Chemicalmolecule` class.

        """
        super().__init__()
        self.chemblid: str = chemblid
        self.mol_data: chembl_webresource_client.query_set.QuerySet = (
            self._get_mol_data())
        self.conformer: Chem.rdchem.Mol = self._get_conformer()
        if self.conformer is not None:
            self.structure: nglview.Structure = nglview.adaptor.RdkitStructure(
                self.conformer)

    def _get_mol_data(self) -> chembl_webresource_client.query_set.QuerySet:
        """ Get the Chemicalmolecule structure from the ChEMBL database
        using the chembl_webresource_client.new_client.molecule API.

        :return: A chembl_webresource_client.query_set.QuerySet containing  
            the 'molecule_chembl_id' and 'molecule_structures' information.

        """
        return new_client.molecule.filter(chembl_id=self.chemblid).only(
            ['molecule_chembl_id', 'molecule_structures'])

    def _get_conformer(self) -> Union[Chem.rdchem.Mol, None]:
        """ Get the 3D conformer molecule structure from the 2D
        structure specified in the 'canonical_smiles' field of self.mol_data.

        :raise: TypeError or KeyError if no SMILES structure is available
            for the specified chemblid.
        :return: A 3D conformer

        """
        try:
            mol: Chem.rdchem.Mol = Chem.AddHs(
                Chem.MolFromSmiles(
                    self.mol_data[0]['molecule_structures']['canonical_smiles']
                ))
            AllChem.EmbedMultipleConfs(mol, useExpTorsionAnglePrefs=True,
                                           useBasicKnowledge=True)
        except (TypeError, KeyError):
            print("No SMILES structure is available for this chemblid.")
            mol = None
        return mol

    def save(self, file_path: Union[str, Path]) -> Path:
        """ Save the ChemicalMolecule 3D structure using the .sdf file format.

        :param file_path: The path to be used when saving the file.
        :raise: FileNotFoundError or OSError if the specified file path
            does not exist
        :return: The path to the saved file.

        """
        if self.conformer is not None:
            try:
                with Chem.SDWriter(str(file_path)) as writer:
                    for cid in range(self.conformer.GetNumConformers()):
                        writer.write(self.conformer, confId=cid)
            except (FileNotFoundError, OSError):
                print("File error: File path does not exist.")
        else:
            print("This ChemicalMolecule cannot be saved as no SMILES "
                  "structure is available for this chemblid.")
        return Path(file_path)
