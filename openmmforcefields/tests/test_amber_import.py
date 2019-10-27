"""
Test AMBER forcefield imports.

"""

import os
import glob
import pytest

from openmmforcefields.utils import get_ffxml_path
amber_ffxml_filenames = glob.glob(os.path.join(get_ffxml_path(), 'amber', '*.xml'))

@pytest.mark.parametrize("filename", amber_ffxml_filenames, ids=lambda filename : f'Importing ffxml file {filename}')
def test_ffxml_import(filename):
    """
    Attempt to load OpenMM ffxml forcefield file.

    Parameters
    ----------
    filename : str
        The AMBER forcefield filename.

    """
    from simtk.openmm import app
    ff = app.ForceField(filename)

def check_ffxml_parameterize(pdb_filename, ffxml_filename):
    """
    Attempt to load OpenMM ffxml forcefield file and parameterize a molecule.

    Parameters
    ----------
    pdb_filename : str
        The PDB filename.
    ffxml_filename : str
        The ffxml forcefield filename.

    """
    from simtk.openmm import app
    pdbfile = app.PDBFile(pdb_filename)
    ff = app.ForceField(ffxml_filename)

def test_amber_import_ff94():
    """
    Test import of ff94

    """
    from simtk.openmm import app
    ff = app.ForceField('amber/ff94.xml')

def test_amber_parameterize_ff94():
    """
    Test parameterizing explicit villin with ff94

    """
    from pkg_resources import resource_filename
    pdb_filename = resource_filename('simtk.openmm.app', 'data/test.pdb')
    check_ffxml_parameterize(pdb_filename, 'amber/ff94.xml')
