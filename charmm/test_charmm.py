from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed import openmm
import glob
import yaml
from collections import OrderedDict
import hashlib
import os
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
import argparse
import csv
import logging
import warnings

def main():
    global verbose
    global no_log
    global logger
    # Set up parser
    parser = argparse.ArgumentParser(description='CHARMM --> OpenMM forcefield conversion test script')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='turns verbosity on')
    parser.add_argument('--no-log', action='store_true',
                        help='turns logging of energies to log.csv off')
    args = parser.parse_args()
    verbose = args.verbose
    no_log = args.no_log

    if not no_log: logger = Logger('log.csv')
    test_charmm()
    if not no_log: logger.close()

def test_charmm():
    """
    Test CHARMM ffxml conversion by computing energy discrepancies between (pdb, psf, toppar) loaded via ParmEd and (pdb, ffxml) loaded via OpenMM ForceField

    """
    # Test systems
    # TODO: Add more test systems generated with CHARMM-GUI.
    testsystems = [
        # name, PDB filename, PSF filename, ffxml filenames, CHARMM toppar filenames
        ('methanol with ions', 'tests/methanol_ions.pdb', 'tests/methanol_ions.psf', ['ffxml/charmm36.xml'], ['toppar/par_all36_cgenff.prm', 'toppar/top_all36_cgenff.rtf','toppar/toppar_water_ions.str']),
    ]

    for (name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames) in testsystems:
        print('Testing %s' % name)
        compare_energies(name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames)

def compare_energies(name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames, system_kwargs=None, tolerance=1e-5, units=u.kilojoules_per_mole):
    """
    Compare energies between (pdb, psf, toppar) loaded via ParmEd and (pdb, ffxml) loaded by OpenMM ForceField

    Parameters
    ----------
    name : str
        Name of the test system
    pdb_filename : str
        Name of PDB file that should contain CRYST entry and PDB format compliant CONECT records for HETATM residues.
    psf_filename : str
        CHARMM PSF file
    ffxml_filenames : list of str
        List of OpenMM ffxml files
    toppar_filenames : list of CHARMM toppar filenames to load into CharmmParameterSet
        List of CHARMM toppar files
    system_kwargs : dict, optional, default=None
        Keyword arguments to pass to CharmmPsfFile.createSystem() and ForceField.CreateSystem() when constructing System objects for energy comparison
    tolerance : float, optional, default=1e-5
        Relative energy discrepancy tolerance
    units : simtk.unit.Unit
        Unit to use for energy comparison

    """
    # Defaults
    if system_kwargs is None:
        system_kwargs = {
            'nonbondedMethod' : app.NoCutoff,
            'constraints' : None,
        }

    # Load PDB file
    pdbfile = app.PDBFile(pdb_filename)

    # Load CHARMM system through ParmEd
    toppar = CharmmParameterSet(*toppar_filenames)
    structure = CharmmPsfFile(psf_filename)
    structure.positions = pdbfile.positions
    system_charmm = structure.createSystem(toppar, **system_kwargs)
    charmm_energies = openmm.energy_decomposition_system(structure, system_charmm, nrg=units)

    # OpenMM system with ffxml
    ff = app.ForceField(*ffxml_filenames)
    system_openmm = ff.createSystem(pdbfile.topology, **system_kwargs)
    topology = openmm.load_topology(pdbfile.topology, system_openmm, xyz=pdbfile.positions)
    omm_energies = openmm.energy_decomposition_system(topology, system_openmm, nrg=units)

    print('charmm_energies')
    print(charmm_energies)
    print('omm_energies')
    print(omm_energies)

    # TODO: Check if discrepancies are larger than tolerance

class Logger():
    # logs testing energies into csv
    def __init__(self, log_file):
        csvfile = open(log_file, 'w')
        fieldnames = ['ffxml_name', 'data_type', 'test_system', 'units', 'HarmonicBondForce',
                      'HarmonicAngleForce', 'PeriodicTorsionForce_dihedrals',
                      'PeriodicTorsionForce_impropers', 'NonbondedForce']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        self.csvfile = csvfile
        self.writer = writer

    def close(self):
        self.csvfile.close()

    def log(self, energies):
        self.writer.writerow(energies)

if __name__ == '__main__':
    main()
