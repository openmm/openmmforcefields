"""Exercise the optional GAFF-only ParmEd dependency."""

import shutil
import subprocess
import sys

import pytest

from openmmforcefields.generators import GAFFTemplateGenerator


def test_generators_import_without_parmed():
    """Normal imports must not require the optional ParmEd dependency for GAFF."""
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; sys.modules['parmed'] = None; "
            "from openmmforcefields.generators import "
            "GAFFTemplateGenerator, SMIRNOFFTemplateGenerator, SystemGenerator",
        ],
        check=True,
        capture_output=True,
        text=True,
    )


def test_gaff_missing_parmed(monkeypatch):
    """Report the missing dependency before starting parameter assignment."""
    monkeypatch.setattr(shutil, "which", lambda _: "/usr/bin/mock-ambertools")
    monkeypatch.setitem(sys.modules, "parmed", None)
    generator = GAFFTemplateGenerator(forcefield="gaff-2.2.20")
    with pytest.raises(RuntimeError, match="ParmEd is required for GAFF"):
        generator.generate_residue_template(object())


def test_gaff_missing_ambertools_still_reports_ambertools(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _: None)
    monkeypatch.setitem(sys.modules, "parmed", None)
    generator = GAFFTemplateGenerator(forcefield="gaff-2.2.20")
    with pytest.raises(RuntimeError, match="AmberTools is required for GAFF"):
        generator.generate_residue_template(object())
