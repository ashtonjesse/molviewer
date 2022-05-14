# Molecular Data Retrieval and Visualisation

## Saving .pdb files
While the PDBList class of the Biopython package can retrieve .pdb files from the www.rcsb.org ftp server, these files are retrieved with the '.ent' extension. Instead, we opted to retrieve .pdb files directly from 'https://files.rcsb.org/download/' as these have the correct '.pdb' extension. Hence the final release of this package does not include Biopython in its dependencies.  
