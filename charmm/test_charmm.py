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

# define NEARLYZERO to replace numerical comparisons to zero
NEARLYZERO = 1e-10

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
        # three-site water models
        ('waterbox TIP3P', 'tests/waterboxes/waterbox-3-site.pdb', 'tests/waterboxes/waterbox-3-site.psf', ['ffxml/waters_ions_default.xml'], ['toppar/toppar_water_ions.str']),
        ('waterbox SPC', 'tests/waterboxes/waterbox-3-site.pdb', 'tests/waterboxes/waterbox-3-site.psf', ['ffxml/waters_ions_spc.xml'], ['toppar/non_charmm/toppar_water_ions_spc.str']),
        ('waterbox SPC/E', 'tests/waterboxes/waterbox-3-site.pdb', 'tests/waterboxes/waterbox-3-site.psf', ['ffxml/waters_ions_spc_e.xml'], ['toppar/non_charmm/toppar_water_ions_spc_e.str']),
        ('waterbox TIP3P PME B', 'tests/waterboxes/waterbox-3-site.pdb', 'tests/waterboxes/waterbox-3-site.psf', ['ffxml/waters_ions_tip3p_pme_b.xml'], ['toppar/non_charmm/toppar_water_ions_tip3p_pme_b.str']),
        ('waterbox TIP3P PME F', 'tests/waterboxes/waterbox-3-site.pdb', 'tests/waterboxes/waterbox-3-site.psf', ['ffxml/waters_ions_tip3p_pme_f.xml'], ['toppar/non_charmm/toppar_water_ions_tip3p_pme_f.str']),
        # multi-site water models
        #('waterbox TIP4P', 'tests/waterboxes/waterbox-4-site.pdb', 'tests/waterboxes/waterbox-4-site.psf', ['ffxml/waters_ions_tip4p.xml'], ['toppar/non_charmm/toppar_water_ions_tip4p.str']),
        #('waterbox TIP4P 2005', 'tests/waterboxes/waterbox-4-site.pdb', 'tests/waterboxes/waterbox-4-site.psf', ['ffxml/waters_ions_tip4p_2005.xml'], ['toppar/non_charmm/toppar_water_ions_tip4p_2005.str']),
        #('waterbox TIP4P-Ew', 'tests/waterboxes/waterbox-4-site.pdb', 'tests/waterboxes/waterbox-4-site.psf', ['ffxml/waters_ions_tip4p_ew.xml'], ['toppar/non_charmm/toppar_water_ions_tip4p_ew.str']),
        #('waterbox TIP5P', 'tests/waterboxes/waterbox-5-site.pdb', 'tests/waterboxes/waterbox-5-site.psf', ['ffxml/waters_ions_tip5p.xml'], ['toppar/non_charmm/toppar_water_ions_tip5p.str']),
        #('waterbox TIP5P-Ew', 'tests/waterboxes/waterbox-5-site.pdb', 'tests/waterboxes/waterbox-5-site.psf', ['ffxml/waters_ions_tip5p_ew.xml'], ['toppar/non_charmm/toppar_water_ions_tip5p_ew.str']),
        # CGenFF
        ('methanol with ions', 'tests/methanol_ions.pdb', 'tests/methanol_ions.psf', ['ffxml/charmm36.xml'], ['toppar/par_all36_cgenff.prm', 'toppar/top_all36_cgenff.rtf','toppar/toppar_water_ions.str']),
    ]

    for (name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames) in testsystems:
        print('')
        print('Testing %s' % name)
        compare_energies(name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames)

def compare_energies(system_name, pdb_filename, psf_filename, ffxml_filenames, toppar_filenames, system_kwargs=None, tolerance=1e-5, units=u.kilojoules_per_mole):
    """
    Compare energies between (pdb, psf, toppar) loaded via ParmEd and (pdb, ffxml) loaded by OpenMM ForceField

    Parameters
    ----------
    system_name : str
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

    # Load CHARMM system through OpenMM
    openmm_toppar = app.CharmmParameterSet(*toppar_filenames)
    openmm_psf = app.CharmmPsfFile(psf_filename)
    openmm_system = openmm_psf.createSystem(openmm_toppar, **system_kwargs)
    with open('system_openmm_charmm_reader.xml', 'w') as f:
        f.write(mm.XmlSerializer.serialize(openmm_system))
    integrator = mm.VerletIntegrator(1.0)
    context = mm.Context(openmm_system, integrator)
    context.setPositions(pdbfile.positions)
    openmm_total_energy = context.getState(getEnergy=True).getPotentialEnergy() / units
    del context, integrator

    # Load CHARMM system through ParmEd
    toppar = CharmmParameterSet(*toppar_filenames)
    structure = CharmmPsfFile(psf_filename)
    #structure.load_parameters(toppar)
    structure.positions = pdbfile.positions
    system_charmm = structure.createSystem(toppar, **system_kwargs)
    with open('system_charmm.xml', 'w') as f:
        f.write(mm.XmlSerializer.serialize(system_charmm))
    charmm_energies = openmm.energy_decomposition_system(structure, system_charmm, nrg=units)
    charmm_total_energy = sum([element[1] for element in charmm_energies])

    # OpenMM system with ffxml
    ff = app.ForceField(*ffxml_filenames)
    system_openmm = ff.createSystem(pdbfile.topology, **system_kwargs)
    print('OpenMM system HarmonicBondForceEnergies')
    with open('system_openmm.xml', 'w') as f:
        f.write(mm.XmlSerializer.serialize(system_openmm))
    topology = openmm.load_topology(pdbfile.topology, system_openmm, xyz=pdbfile.positions)
    omm_energies = openmm.energy_decomposition_system(topology, system_openmm, nrg=units)
    ffxml_tootal_energy = sum([element[1] for element in omm_energies])

    print('OpenMM CHARMM reader total energy: %f' % openmm_total_energy)
    print('ParmEd CHARMM reader total energy: %f' % charmm_total_energy)
    print('OPENMM ffxml total energy: %f' % openmm_total_energy)
    print('TOTAL ERROR: %f' % (openmm_total_energy - charmm_total_energy))

    print('ParmEd CHARMM reader energy decomposition:')
    print(charmm_energies)
    print('OpenMM ffxml ForceField energy decomposition:')
    print(omm_energies)

    # TODO : Automate comparison
    return

    # calc rel energies and assert
    rel_energies = []
    for i, j in zip(charmm_energies, omm_energies):
        if i[0] != j[0]:
            raise Exception('Mismatch in energy tuples naming.')
        if abs(i[1]) > NEARLYZERO:
            rel_energies.append((i[0], abs((i[1]-j[1])/i[1])))
        else:
            if abs(j[1]) > NEARLYZERO:
                raise AssertionError('One of the CHARMM %s energies (%s) for %s is zero, '
                      'while the corresponding OpenMM energy is non-zero' %
                      (system_name, i[0], ffxml))
            rel_energies.append((i[0], 0))

    dihedrals_done = False
    for i in rel_energies:
        if i[0] != 'PeriodicTorsionForce':
            if i[1] > tolerance:
                raise AssertionError('%s energies (%s, %f) outside of allowed tolerance (%f) for %s' %
                                     (system_name, i[0], i[1], tolerance, ffxml))
        else:
            if not dihedrals_done:
                if i[1] > tolerance:
                    raise AssertionError('%s energies (%s, %f) outside of allowed tolerance (%f) for %s' %
                                         (system_name, i[0], i[1], tolerance, ffxml))
                dihedrals_done = True
            else: #impropers
                if i[1] > improper_tolerance:
                    raise AssertionError('%s energies (%s-impropers, %f) outside of allowed tolerance (%f) for %s' %
                                         (system_name, i[0], i[1], improper_tolerance, ffxml))

    # logging
    if not no_log:
        charmm_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        charmm_energies_log['ffxml_name'] = ffxml
        charmm_energies_log['test_system'] = system_name
        charmm_energies_log['data_type'] = 'CHARMM'
        charmm_energies_log['units'] = units
        omm_energies_log['ffxml_name'] = ffxml
        omm_energies_log['test_system'] = system_name
        omm_energies_log['data_type'] = 'OpenMM'
        omm_energies_log['units'] = units
        rel_energies_log['ffxml_name'] = ffxml
        rel_energies_log['test_system'] = system_name
        rel_energies_log['data_type'] = 'abs((CHARMM-OpenMM)/CHARMM)'
        dihedrals_done = False
        for item in amber_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                charmm_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                charmm_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                charmm_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in omm_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                omm_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in rel_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                rel_energies_log[item[0]] = item[1]

        logger.log(charmm_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)

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
