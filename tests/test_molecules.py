from pathlib import Path
from rdkit import Chem
from src.molviewer.molecules import ChemicalMolecule, Macromolecule
import pytest
import nglview


@pytest.fixture()
def protease():
    return Macromolecule(pdbid="6LU7")


def test_setting_pdbid(protease):
    assert protease.pdbid == "6LU7"


def test_show_macromol_existing_viewer(protease):
    viewer = protease.show(nglview.NGLWidget())
    assert isinstance(viewer, nglview.NGLWidget)


def test_show_macromol_no_existing_viewer(protease):
    viewer = protease.show()
    assert isinstance(viewer, nglview.NGLWidget)


def test_show_macromol_wrong_existing_viewer(protease):
    viewer = protease.show(
        'an argument which is not of type nglview.NGLWidget')
    assert isinstance(viewer, str)


def test_save_macromol_correct_file_path(protease):
    file_path = protease.save(Path.cwd().joinpath("tests","6LU7.pdb"))
    assert file_path.is_file()


def test_save_macromol_incorrect_file_path(protease):
    file_path = protease.save(Path.cwd().joinpath("nofolder","6LU7.pdb"))
    assert not file_path.is_file()


@pytest.fixture()
def lopinavir():
    return ChemicalMolecule(chemblid="CHEMBL729")


def test_setting_chemblid(lopinavir):
    assert lopinavir.chemblid == "CHEMBL729"


def test_retrieve_mol_data(lopinavir):
    assert len(lopinavir.mol_data) == 1
    assert list(lopinavir.mol_data[0].keys()) == [
        'molecule_chembl_id', 'molecule_structures']
    assert lopinavir.mol_data[0]['molecule_chembl_id'] == "CHEMBL729"
    assert 'canonical_smiles' in lopinavir.mol_data[0][
        'molecule_structures'].keys()


def test_conformer_generation(lopinavir):
    assert isinstance(lopinavir.conformer, Chem.rdchem.Mol)


def test_show_chemmol_existing_viewer(lopinavir):
    viewer = lopinavir.show(nglview.NGLWidget())
    assert isinstance(viewer, nglview.NGLWidget)


def test_show_chemmol_no_existing_viewer(lopinavir):
    viewer = lopinavir.show()
    assert isinstance(viewer, nglview.NGLWidget)


def test_show_chemmol_wrong_existing_viewer(lopinavir):
    viewer = lopinavir.show(
        'an argument which is not of type nglview.NGLWidget')
    assert isinstance(viewer, str)


def test_save_chemmol_correct_file_path(lopinavir):
    file_path = lopinavir.save(Path.cwd().joinpath("tests","CHEMBL729.sdf"))
    assert file_path.is_file()


def test_save_chemmol_incorrect_file_path(lopinavir):
    file_path = lopinavir.save(Path.cwd().joinpath("nofolder","CHEMBL729.sdf"))
    assert not file_path.is_file()

@pytest.fixture()
def erbB():
    return ChemicalMolecule(chemblid="CHEMBL3848")


def test_chemblid_with_no_smiles(erbB):
    assert erbB.conformer is None


def test_show_chemblid_with_no_smiles(erbB):
    viewer = erbB.show()
    assert isinstance(viewer, nglview.NGLWidget)


def test_save_chemblid_with_no_smiles(erbB):
    file_path = erbB.save(Path.cwd().joinpath("tests","CHEMBL3848.sdf"))
    assert not file_path.is_file()
