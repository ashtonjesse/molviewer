Introduction
============

``molviewer`` is a Python package which aims to provide an easy and
intuitive way of retrieving and viewing 3D models of molecules from
the `rcsb <https://www.rcsb.org>`_ and `ChEMBL <https://www.ebi.ac
.uk/chembl/>`_ databases. It is intended to be used within a
``Jupyter notebook`` and essentially provides an extension
of the ``nglview`` Python package (see `here <https://github
.com/nglviewer/nglview>`_).

The aim here was to define two classes: :class:`Macromolecule` for molecules
from the rcsb database, and :class:`ChemicalMolecule` for small molecules from the
ChEMBL database. Both classes wrap the ``nglview.NGLWidget``
methods for viewing multiple 3D molecule models and enable the user to save
the molecular information to file. In the case of :class:`ChemicalMolecule`, 2D
SMILES structure is converted in to a 3D conformer using methods
available in the `rdkit <https://www.rdkit.org/>`_ Python package.

The package is built using Python 3.8+ and implements typing
conventions according to PEP484 and PEP526.

Motivation
**********

By providing methods for retrieving and visualizing the 3D structure of
molecules, we enable detailed examination of how molecules interact.
