"""
Test AMBER forcefield imports.

"""

import os
import glob
from functools import partial

def check_ffxml_import(filename):
    """
    Attempt to load OpenMM ffxml forcefield file.

    Parameters
    ----------
    filename : str
        The AMBER forcefield filename.

    """
    from simtk.openmm import app
    ff = app.ForceField(filename)

def test_amber_import_ff94():
    """
    Test import of ff94

    """
    from simtk.openmm import app
    ff = app.ForceField('amber/ff94.xml')

def test_all_amber_imports():
    """
    Test all available AMBER ffxml files.

    """
    from openmmforcefields.utils import get_ffxml_path
    amber_ffxml_filenames = glob.glob(os.path.join(get_ffxml_path(), 'amber', '*.xml'))
    # Try to import all ffxml files
    for filename in amber_ffxml_filenames:
        f = partial(check_ffxml_import, filename)
        f.description = "Importing ffxml file '%s'" % filename
        yield f
