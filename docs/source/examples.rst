Examples
=============

Installation/Usage:
*******************
Create a virtual environment with Python 3.8+ interpreter and activate it.
For Windows users this can be achieved with the following commands:

.. code-block:: python

    -m venv \path\to\new\virtual\environment
    .\Scripts\activate

If loading scripts is disabled then first run this command in the Windows
powershell:

.. code-block:: python

    Set-ExecutionPolicy Unrestricted -Scope Process

Then install the ``molviewer`` package:

.. code-block:: python

    py -m pip install -i https://test.pypi.org/pypi/ --extra-index-url https://pypi.org/simple molviewer==0.0.13

Run ``Jupyter`` with the following command:


.. code-block:: python

    jupyter notebook

Create a new ``notebook`` and import the ``molviewer`` package:

.. code-block:: python

    from molviewer.molecules import ChemicalMolecule, Macromolecule

Macromolecules can then be loaded from the `rcsb <https://www.rcsb.org>`_
database by specifying a valid PDB ID code: a 4-character alphanumeric
string, e.g.:

.. code-block:: python

    protease = Macromolecule(pdbid="6LU7")

View the 3D molecule model with the following commands:

.. code-block:: python

    viewer = protease.show()
    viewer

A small chemical molecule can then be loaded from the `ChEMBL
<https://www.ebi.ac.uk/chembl/>`_ database by specifying a valid CheMBL ID
code (for a molecule with SMILES structure in the database), e.g.:

.. code-block:: python

    lopinavir = ChemicalMolecule(chemblid="CHEMBL729")

The 3D model can then be added to the current viewer with the following
command:

.. code-block:: python

    viewer = lopinavir.show(viewer)

Both molecules can then be saved with the following commands:

.. code-block:: python

    lopinavir.save('valid\\path\\to\\file')
    protease.save('valid\\path\\to\\file')
