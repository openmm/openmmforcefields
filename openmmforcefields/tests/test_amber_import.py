"""
Test AMBER forcefield imports.

"""

import pathlib

import pytest

from openmmforcefields.utils import get_ffxml_path

amber_ffxml_filenames: list[str] = ["amber/" + file.name for file in pathlib.Path(get_ffxml_path()).glob("amber/*xml")]


@pytest.mark.parametrize(
    "filename",
    amber_ffxml_filenames,
    ids=lambda filename: f"Importing ffxml file {filename}",
)
def test_ffxml_import(filename):
    """
    Attempt to load OpenMM ffxml forcefield file.

    Parameters
    ----------
    filename : str
        The AMBER forcefield filename.

    """
    from openmm import app

    # Handle special cases
    if filename == "amber/phosaa10.xml":
        # Must be used with ff99SB.xml
        filenames = ["amber/ff99SB.xml", "amber/phosaa10.xml"]
        app.ForceField(*filenames)
    elif filename == "amber/phosaa14SB.xml":
        # Must be used with ff14SB.xml
        filenames = ["amber/ff14SB.xml", "amber/phosaa14SB.xml"]
        app.ForceField(*filenames)
    elif filename == "amber/GLYCAM_06j-1.xml":
        # Must be used with protein.ff14SB.xml
        filenames = ["amber/protein.ff14SB.xml", "amber/GLYCAM_06j-1.xml"]
        app.ForceField(*filenames)
    else:
        app.ForceField(filename)


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
    from openmm import app

    app.PDBFile(pdb_filename)
    app.ForceField(ffxml_filename)


def test_amber_import_ff94():
    """
    Test import of ff94

    """
    from openmm import app

    app.ForceField("amber/ff94.xml")


def test_amber_parameterize_ff94():
    """
    Test parameterizing explicit villin with ff94

    """
    import importlib_resources

    pdb_filename = str(importlib_resources.files("openmm.app") / "data" / "test.pdb")
    check_ffxml_parameterize(pdb_filename, "amber/ff94.xml")
