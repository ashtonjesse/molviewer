from molecules import ChemicalMolecule, Macromolecule # type: ignore

protease = Macromolecule(pdbid="6LU7")
lopinavir = ChemicalMolecule(chemblid="CHEMBL729")
viewer = protease.show()
lopinavir.show(viewer)
file_path = protease.save("C:\\Users\\jash042\\Downloads\\6LU7.pdb")
file_path = lopinavir.save("C:\\Users\\jash042\\Downloa\\CHEMBL729.sdf")
my_test = 1
