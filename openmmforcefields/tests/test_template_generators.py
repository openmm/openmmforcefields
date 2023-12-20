import copy
import logging
import os
import pytest
import tempfile
import unittest
import pytest

import numpy as np
import openmm
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField as OFFForceField
from openff.units import unit as OFFUnit
from openmm.app import PME, ForceField, Modeller, NoCutoff, PDBFile

from openmmforcefields.generators import (
    EspalomaTemplateGenerator,
    GAFFTemplateGenerator,
    SMIRNOFFTemplateGenerator,
)
from openmmforcefields.utils import get_data_filename

_logger = logging.getLogger("openmmforcefields.tests.test_template_generators")

CI = ('CI' in os.environ)

################################################################################
# Tests
################################################################################

class TestGAFFTemplateGenerator(unittest.TestCase):
    TEMPLATE_GENERATOR = GAFFTemplateGenerator

    amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']

    def filter_molecules(self, molecules):
        """
        Filter molecules to speed up tests, especially on travis.

        Parameters
        ----------
        molecules : list of openff.toolkit.topology.Molecule
            The input list of molecules to be filtered

        Returns
        -------
        molecules : list of openff.toolkit.topology.Molecule
            The filtered list of molecules to be filtered

        """
        # TODO: Eliminate molecules without fully-specified stereochemistry
        # Select some small molecules for fast testing
        MAX_ATOMS = 40
        molecules = [ molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS ]
        # Cut down number of tests for continuous integration
        MAX_MOLECULES = 50 if not CI else 4
        molecules = molecules[:MAX_MOLECULES]

        return molecules

    def setUp(self):
        from openff.units import unit
        from openff.units.openmm import ensure_quantity

        # TODO: Harmonize with test_system_generator.py infrastructure
        # Read test molecules
        filename = get_data_filename("minidrugbank/MiniDrugBank-without-unspecified-stereochemistry.sdf")
        molecules = Molecule.from_file(filename, allow_undefined_stereo=True)

        # DEBUG: Insert acetone perturbed from planarity as first test molecule, since it fails quickly if something is wrong
        molecule = Molecule.from_smiles('C=O')
        molecule.generate_conformers(n_conformers=1)

        uses_old_api = hasattr(molecule.atoms[0], "element")

        molecule.conformers[0][0,0] += ensure_quantity(
            unit.Quantity(0.1, unit.angstroms),
            "openmm" if uses_old_api else "openff",
        )

        molecules.insert(0, molecule)
        # DEBUG END

        # Filter molecules as appropriate
        self.molecules = self.filter_molecules(molecules)

        # Suppress DEBUG logging from various packages
        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

    def test_version(self):
        """Test version"""
        for forcefield in GAFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator = GAFFTemplateGenerator(forcefield=forcefield)
            import re
            result = re.match(r'^gaff-(?P<major_version>\d+)\.(?P<minor_version>\d+)$', forcefield)
            assert generator.forcefield == forcefield
            assert generator.gaff_version == result['major_version'] + '.' + result['minor_version']
            assert generator.gaff_major_version == result['major_version']
            assert generator.gaff_minor_version == result['minor_version']
            assert generator.gaff_dat_filename.endswith(forcefield + '.dat')
            assert os.path.exists(generator.gaff_dat_filename)
            assert generator.gaff_xml_filename.endswith(forcefield + '.xml')
            assert os.path.exists(generator.gaff_xml_filename)

    def test_create(self):
        """Test template generator creation"""
        # Create an empty generator
        generator = self.TEMPLATE_GENERATOR()
        # Create a generator that knows about a few molecules
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a new database file
            generator = self.TEMPLATE_GENERATOR(molecules=self.molecules, cache=cache)
            del generator
            # Reopen it (with cache still empty)
            generator = self.TEMPLATE_GENERATOR(molecules=self.molecules, cache=cache)
            del generator

    def test_add_molecules(self):
        """Test that molecules can be added to template generator after its creation"""
        # Create a generator that does not know about any molecules
        generator = self.TEMPLATE_GENERATOR()
        # Create a ForceField
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)

        # Check that parameterizing a molecule fails
        molecule = self.molecules[0]
        try:
            # This should fail with an exception
            openmm_topology = molecule.to_topology().to_openmm()
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
        except ValueError as e:
            # Exception 'No template found...' is expected
            assert str(e).startswith('No template found')

        # Now add the molecule to the generator and ensure parameterization passes
        generator.add_molecules(molecule)
        openmm_topology = molecule.to_topology().to_openmm()
        try:
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
        except Exception as e:
            print(forcefield._atomTypes.keys())
            from openff.units.openmm import ensure_quantity
            PDBFile.writeFile(
                openmm_topology,
                ensure_quantity(molecule.conformers[0], "openmm"),
            )
            raise e
        assert system.getNumParticles() == molecule.n_atoms

        # Add multiple molecules, including repeats
        generator.add_molecules(self.molecules)

        # Ensure all molecules can be parameterized
        for molecule in self.molecules:
            openmm_topology = molecule.to_topology().to_openmm()
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert system.getNumParticles() == molecule.n_atoms

    def charges_from_system(self, system):
        """Extract dimensionless partial charges from a System

        Parameters
        ----------
        system : openmm.System
            The System from which partial charges are to be extracted

        Returns
        -------
        charges : np.array of shape [n_atoms]
            The dimensionless partial charges (implicitly in units of elementary_charge)

        """
        from openmm import unit
        system_charges = list()
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        for particle_index in range(system.getNumParticles()):
            charge, sigma, epsilon = forces['NonbondedForce'].getParticleParameters(particle_index)
            system_charges.append(charge / unit.elementary_charge)
        system_charges = np.array(system_charges)

        return system_charges

    def charges_are_equal(self, system, molecule):
        """Return True if partial charges are equal

        Parameters
        ----------
        system : openmm.System
            The System from which partial charges are to be extracted
        molecule : openmmforcefield.topology.Molecule
            The Molecule from which partial charges are to be extracted

        Returns
        -------
        result : bool
            True if the partial charges are equal, False if not
        """
        from openff.units import unit
        from openff.units.openmm import ensure_quantity

        assert system.getNumParticles() == molecule.n_atoms

        # charges_from_system returns a NumPy array that we trust to be implicitly e
        system_charges: np.ndarray = self.charges_from_system(system)

        # type(molecule.partial_charges) depends on the toolkit version
        molecule_charges: np.ndarray = ensure_quantity(
            molecule.partial_charges, "openff",
        ).m_as(unit.elementary_charge)

        result = np.allclose(system_charges, molecule_charges)

        if not result:
            _logger.info('Charges are not equal')
            _logger.info(f'system charges  : {system_charges}')
            _logger.info(f'molecule charges: {molecule_charges}')

        return result

    def test_charge(self):
        """Test that charges are nonzero after charging if the molecule does not contain user charges"""
        # Create a generator that does not know about any molecules
        generator = self.TEMPLATE_GENERATOR()
        # Create a ForceField
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)

        # Check that parameterizing a molecule using user-provided charges produces expected charges
        from openmm import unit
        molecule = self.molecules[0]
        # Ensure partial charges are initially zero
        assert (molecule.partial_charges is None) or np.all(molecule.partial_charges / unit.elementary_charge == 0)
        # Add the molecule
        generator.add_molecules(molecule)
        # Create the System
        openmm_topology = molecule.to_topology().to_openmm()
        system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
        # Ensure charges are no longer zero
        assert not np.all(self.charges_from_system(system) == 0), "System has zero charges despite molecule not being charged"

    def test_charge_from_molecules(self):
        """Test that user-specified partial charges are used if requested"""
        from openff.units.openmm import ensure_quantity

        # Create a generator that does not know about any molecules
        generator = self.TEMPLATE_GENERATOR()
        # Create a ForceField
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)

        # Check that parameterizing a molecule using user-provided charges produces expected charges

        molecule = self.molecules[0]
        uses_old_api = hasattr(molecule.atoms[0], "element")

        if uses_old_api:
            from openmm import unit
        else:
            from openff.units import unit

        # Populate the molecule with arbitrary partial charges that still sum to 0.0
        molecule.partial_charges = unit.Quantity(
            np.linspace(-0.5, 0.5, molecule.n_atoms),
            unit.elementary_charge,
        )

        assert (molecule.partial_charges is not None)

        assert not np.all(ensure_quantity(molecule.partial_charges, "openff").m == 0)

        generator.add_molecules(molecule)

        system = forcefield.createSystem(molecule.to_topology().to_openmm(), nonbondedMethod=NoCutoff)

        assert self.charges_are_equal(system, molecule)

    def test_debug_ffxml(self):
        """Test that debug ffxml file is created when requested"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            debug_ffxml_filename = os.path.join(tmpdirname, 'molecule.ffxml')
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a generator that only knows about one molecule
            molecule = self.molecules[0]
            generator = self.TEMPLATE_GENERATOR(molecules=molecule, cache=cache)
            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Ensure no file is created
            openmm_topology = molecule.to_topology().to_openmm()
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert not os.path.exists(debug_ffxml_filename)
            # Enable debug file output creation
            forcefield = ForceField()
            forcefield.registerTemplateGenerator(generator.generator)
            generator.debug_ffxml_filename = debug_ffxml_filename
            # Ensure that an ffxml file is created
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert os.path.exists(debug_ffxml_filename)
            # Ensure we can use that file to create a new force field
            forcefield_from_ffxml = ForceField()
            if hasattr(generator, 'gaff_xml_filename'):
                forcefield_from_ffxml.loadFile(generator.gaff_xml_filename)
            forcefield_from_ffxml.loadFile(debug_ffxml_filename)
            system2 = forcefield_from_ffxml.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            # TODO: Test that systems are equivalent
            assert system.getNumParticles() == system2.getNumParticles()

    def test_cache(self):
        """Test template generator cache capability"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a generator that also has a database cache
            cache = os.path.join(tmpdirname, 'db.json')
            generator = self.TEMPLATE_GENERATOR(molecules=self.molecules, cache=cache)
            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules
            for molecule in self.molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)

            # Check database contents
            def check_cache(generator, n_expected):
                """
                Check database contains number of expected records

                Parameters
                ----------
                generator : SmallMoleculeTemplateGenerator
                    The generator whose cache should be examined
                n_expected : int
                    Number of expected records
                """
                from tinydb import TinyDB
                db = TinyDB(generator._cache)
                table = db.table(generator._database_table_name)
                db_entries = table.all()
                db.close()
                n_entries = len(db_entries)
                assert (n_entries == n_expected), \
                    f"Expected {n_expected} entries but database has {n_entries}\n db contents: {db_entries}"

            check_cache(generator, len(self.molecules))

            # Clean up, forcing closure of database
            del forcefield, generator

            # Create a generator that also uses the database cache but has no molecules
            print('Creating new generator with just cache...')
            generator = self.TEMPLATE_GENERATOR(cache=cache)
            # Check database still contains the molecules we expect
            check_cache(generator, len(self.molecules))
            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules; this should succeed
            for molecule in self.molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)

    def test_add_solvent(self):
        """Test using openmm.app.Modeller to add solvent to a small molecule parameterized by template generator"""
        # Select a molecule to add solvent around
        from openff.units.openmm import ensure_quantity
        from openmm import unit
        molecule = self.molecules[0]
        openmm_topology = molecule.to_topology().to_openmm()
        openmm_positions = ensure_quantity(molecule.conformers[0], "openmm")
        # Try adding solvent without residue template generator; this will fail
        forcefield = ForceField('tip3p.xml')
        # Add solvent to a system containing a small molecule
        modeller = Modeller(openmm_topology, openmm_positions)
        try:
            modeller.addSolvent(forcefield, model='tip3p', padding=6.0*unit.angstroms)
        except ValueError as e:
            pass

        # Create a generator that knows about a few molecules
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        # Add to the forcefield object
        forcefield.registerTemplateGenerator(generator.generator)
        # Add solvent to a system containing a small molecule
        # This should succeed
        modeller.addSolvent(forcefield, model='tip3p', padding=6.0*unit.angstroms)

    def test_jacs_ligands(self):
        """Use template generator to parameterize the Schrodinger JACS set of ligands"""
        jacs_systems = {
            #'bace'     : { 'prefix' : 'Bace' },
            #'cdk2'     : { 'prefix' : 'CDK2' },
            'jnk1'     : { 'prefix' : 'Jnk1' },
            'mcl1'     : { 'prefix' : 'MCL1' },
            #'p38'      : { 'prefix' : 'p38' },
            'ptp1b'    : { 'prefix' : 'PTP1B' },
            'thrombin' : { 'prefix' : 'Thrombin' },
            #'tyk2'     : { 'prefix' : 'Tyk2' },
        }
        for system_name in jacs_systems:
            prefix = jacs_systems[system_name]['prefix']
            # Load molecules
            ligand_sdf_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_ligands.sdf'))
            print(f'Reading molecules from {ligand_sdf_filename} ...')
            molecules = Molecule.from_file(ligand_sdf_filename, allow_undefined_stereo=True)
            # Ensure this is a list
            try:
                nmolecules = len(molecules)
            except TypeError:
                molecules = [molecules]

            print(f'Read {len(molecules)} molecules from {ligand_sdf_filename}')
            #molecules = self.filter_molecules(molecules)
            MAX_MOLECULES = len(molecules) if not CI else 3
            molecules = molecules[:MAX_MOLECULES]
            print(f'{len(molecules)} molecules remain after filtering')

            # Create template generator with local cache
            cache = os.path.join(get_data_filename(os.path.join('perses_jacs_systems', system_name)), 'cache.json')
            generator = self.TEMPLATE_GENERATOR(molecules=molecules, cache=cache)

            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)

            # Parameterize all molecules
            print(f'Caching all molecules for {system_name} at {cache} ...')
            n_success = 0
            n_failure = 0
            for molecule in molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                try:
                    forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
                    n_success += 1
                except Exception as e:
                    n_failure += 1
                    print(e)
            print(f'{n_failure}/{n_success+n_failure} ligands failed to parameterize for {system_name}')

    def test_jacs_complexes(self):
        """Use template generator to parameterize the Schrodinger JACS set of complexes"""
        # TODO: Uncomment working systems when we have cleaned up the input files
        jacs_systems = {
            #'bace'     : { 'prefix' : 'Bace' },
            #'cdk2'     : { 'prefix' : 'CDK2' },
            #'jnk1'     : { 'prefix' : 'Jnk1' },
            'mcl1'     : { 'prefix' : 'MCL1' },
            #'p38'      : { 'prefix' : 'p38' },
            #'ptp1b'    : { 'prefix' : 'PTP1B' },
            #'thrombin' : { 'prefix' : 'Thrombin' },
            #'tyk2'     : { 'prefix' : 'Tyk2' },
        }
        for system_name in jacs_systems:
            prefix = jacs_systems[system_name]['prefix']
            # Read molecules
            ligand_sdf_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_ligands.sdf'))
            print(f'Reading molecules from {ligand_sdf_filename} ...')
            molecules = Molecule.from_file(ligand_sdf_filename, allow_undefined_stereo=True)
            try:
                nmolecules = len(molecules)
            except TypeError:
                molecules = [molecules]
            print(f'Read {len(molecules)} molecules from {ligand_sdf_filename}')

            # Read ParmEd Structures
            import parmed
            from openmm import unit
            protein_pdb_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_protein.pdb'))
            print(f'Reading protein from {protein_pdb_filename} ...')
            #protein_structure = parmed.load_file(protein_pdb_filename) # NOTE: This mis-interprets distorted geometry and sequentially-numbered residues that span chain breaks
            pdbfile = PDBFile(protein_pdb_filename)
            protein_structure = parmed.openmm.load_topology(pdbfile.topology, xyz=pdbfile.positions.value_in_unit(unit.angstroms))
            ligand_structures = parmed.load_file(ligand_sdf_filename)
            try:
                nmolecules = len(ligand_structures)
            except TypeError:
                ligand_structures = [ligand_structures]
            assert len(ligand_structures) == len(molecules)

            # Filter molecules
            MAX_MOLECULES = 6 if not CI else 3
            molecules = molecules[:MAX_MOLECULES]
            ligand_structures = ligand_structures[:MAX_MOLECULES]
            print(f'{len(molecules)} molecules remain after filtering')

            # Create complexes
            complex_structures = [ (protein_structure + ligand_structure) for ligand_structure in ligand_structures ]

            # Create template generator with local cache
            cache = os.path.join(get_data_filename(os.path.join('perses_jacs_systems', system_name)), 'cache.json')
            generator = self.TEMPLATE_GENERATOR(molecules=molecules, cache=cache)

            # Create a ForceField
            forcefield = ForceField(*self.amber_forcefields)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)

            # Parameterize all complexes
            print(f'Caching all molecules for {system_name} at {cache} ...')
            for ligand_index, complex_structure in enumerate(complex_structures):
                openmm_topology = complex_structure.topology
                molecule = molecules[ligand_index]

                # Delete hydrogens from terminal protein residues
                # TODO: Fix the input files so we don't need to do this
                from openmm import app
                modeller = app.Modeller(complex_structure.topology, complex_structure.positions)
                residues = [residue for residue in modeller.topology.residues() if residue.name != 'UNL']
                termini_ids = [residues[0].id, residues[-1].id]
                #hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol in ['H'] and atom.residue.name != 'UNL']
                hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol in ['H'] and atom.residue.id in termini_ids]
                modeller.delete(hs)
                modeller.addHydrogens(forcefield)

                # Parameterize protein:ligand complex in vacuum
                print(f' Parameterizing {system_name} : {molecule.to_smiles()} in vacuum...')
                forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

                # Parameterize protein:ligand complex in solvent
                print(f' Parameterizing {system_name} : {molecule.to_smiles()} in explicit solvent...')
                modeller.addSolvent(forcefield, padding=0*unit.angstroms, ionicStrength=300*unit.millimolar)
                forcefield.createSystem(modeller.topology, nonbondedMethod=PME)

    def test_parameterize(self):
        """Test parameterizing molecules with template generator for all supported force fields"""
        # Test all supported small molecule force fields
        for small_molecule_forcefield in self.TEMPLATE_GENERATOR.INSTALLED_FORCEFIELDS:
            if "ff14sb" in small_molecule_forcefield:
                continue
            if "tip" in small_molecule_forcefield:
                continue
            if "opc" in small_molecule_forcefield:
                continue

            print(f'Testing {small_molecule_forcefield}')
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = self.TEMPLATE_GENERATOR(molecules=self.molecules, forcefield=small_molecule_forcefield)
            # Check that we have loaded the right force field
            assert generator.forcefield == small_molecule_forcefield
            # Create a ForceField with the appropriate small molecule force field
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules

            from openmmforcefields.utils import Timer
            for molecule in self.molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                with Timer() as t1:
                    system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == molecule.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == molecule.n_atoms
                assert (t2.interval() < t1.interval())

    def test_multiple_registration(self):
        """Test registering the template generator with multiple force fields"""
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        NUM_FORCEFIELDS = 2 # number of force fields to test
        forcefields = list()
        for index in range(NUM_FORCEFIELDS):
            forcefield = ForceField()
            forcefield.registerTemplateGenerator(generator.generator)
            forcefields.append(forcefield)

        # Parameterize a molecule in each force field instance
        molecule = self.molecules[0]
        openmm_topology = molecule.to_topology().to_openmm()
        for forcefield in forcefields:
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert system.getNumParticles() == molecule.n_atoms

    @staticmethod
    def compute_energy(system, positions):
        """Compute potential energy and Force components for an OpenMM system.

        Parameters
        ----------
        system : openmm.System
            The System object
        positions : openmm.unit.Quantity of shape (nparticles,3) with units compatible with nanometers
            The positions for which energy is to be computed

        Returns
        -------
        openmm_energy : dict of str : openmm.unit.Quantity
            openmm_energy['total'] is the total potential energy
            openmm_energy['components'][forcename] is the potential energy for the specified component force
        openmm_forces : dict of str : openmm.unit.Quantity
            openmm_forces['total'] is the total force
            openmm_forces['components'][forcename] is the force for the specified component force
        """
        system = copy.deepcopy(system)
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
        platform = openmm.Platform.getPlatformByName('Reference')
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        openmm_energy = {
            'total' : context.getState(getEnergy=True).getPotentialEnergy(),
            'components' : { system.getForce(index).__class__.__name__ : context.getState(getEnergy=True, groups=(1 << index)).getPotentialEnergy() for index in range(system.getNumForces()) },
            }

        openmm_forces = {
            'total' : context.getState(getForces=True).getForces(asNumpy=True),
            'components' : { system.getForce(index).__class__.__name__ : context.getState(getForces=True, groups=(1 << index)).getForces(asNumpy=True) for index in range(system.getNumForces()) },
            }

        del context, integrator
        return openmm_energy, openmm_forces

    @classmethod
    def compare_energies(cls, molecule, template_generated_system, reference_system):
        """Compare energies between OpenMM System generated by reference method and OpenMM System generated by ForceField template.

        The OpenMM System object created internally by the reference method is used to
        avoid any issues with stochasticity of partial charges due to conformer generation.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The Molecule object to compare energy components (including positions)
        template_generated_system : openmm.System
            System generated by OpenMM ForceField template
        reference_system : openmm.System
            System generated by reference parmaeterization engine

        """
        from openff.units.openmm import ensure_quantity

        # Compute energies
        reference_energy, reference_forces = cls.compute_energy(
            template_generated_system,
            ensure_quantity(molecule.conformers[0], "openmm"),
        )
        template_energy, template_forces = cls.compute_energy(
            reference_system,
            ensure_quantity(molecule.conformers[0], "openmm"),
        )

        from openmm import unit

        def write_xml(filename, system):
            with open(filename, 'w') as outfile:
                print(f'Writing {filename}...')
                outfile.write(openmm.XmlSerializer.serialize(system))
                # DEBUG
                print(openmm.XmlSerializer.serialize(system))

        # Make sure both systems contain the same energy components
        reference_components = set(reference_energy['components'])
        template_components = set(template_energy['components'])
        if len(reference_components.difference(template_components)) > 0:
            raise Exception(f'Reference system contains components {reference_components.difference(template_components)} that do not appear in template-generated system.')
        if len(template_components.difference(reference_components)) > 0:
            raise Exception(f'Template-generated system contains components {template_components.difference(reference_components)} that do not appear in reference system.')
        components = reference_components

        # Compare energies
        ENERGY_DEVIATION_TOLERANCE = 1.0e-2 * unit.kilocalories_per_mole
        delta = (template_energy['total'] - reference_energy['total'])
        if abs(delta) > ENERGY_DEVIATION_TOLERANCE:
            # Show breakdown by components
            print('Energy components:')
            print(f"{'component':24} {'Template (kcal/mol)':>20} {'Reference (kcal/mol)':>20}")
            for key in components:
                reference_component_energy = reference_energy['components'][key]
                template_component_energy = template_energy['components'][key]
                print(f'{key:24} {(template_component_energy/unit.kilocalories_per_mole):20.3f} {(reference_component_energy/unit.kilocalories_per_mole):20.3f} kcal/mol')
            print(f'{"TOTAL":24} {(template_energy["total"]/unit.kilocalories_per_mole):20.3f} {(reference_energy["total"]/unit.kilocalories_per_mole):20.3f} kcal/mol')
            write_xml('reference_system.xml', reference_system)
            write_xml('template_system.xml', template_system)  # What's this? This variable does not exist
            raise Exception(f'Energy deviation for {molecule.to_smiles()} ({delta/unit.kilocalories_per_mole} kcal/mol) exceeds threshold ({ENERGY_DEVIATION_TOLERANCE})')

        # Compare forces
        def norm(x):
            N = x.shape[0]
            return np.sqrt((1.0/N) * (x**2).sum())
        def relative_deviation(x, y):
            FORCE_UNIT = unit.kilocalories_per_mole / unit.angstroms
            if hasattr(x, 'value_in_unit'):
                x = x / FORCE_UNIT
            if hasattr(y, 'value_in_unit'):
                y = y / FORCE_UNIT

            if norm(y) > 0:
                return norm(x-y) / np.sqrt(norm(x)**2 + norm(y)**2)
            else:
                return 0

        RELATIVE_FORCE_DEVIATION_TOLERANCE = 1.0e-5
        relative_force_deviation = relative_deviation(template_forces['total'], reference_forces['total'])
        if relative_force_deviation > RELATIVE_FORCE_DEVIATION_TOLERANCE:
            # Show breakdown by components
            print('Force components:')
            print(f"{'component':24} {'relative deviation':>24}")
            for key in components:
                print(f"{key:24} {relative_deviation(template_forces['components'][key], reference_forces['components'][key]):24.10f}")
            print(f'{"TOTAL":24} {relative_force_deviation:24.10f}')
            write_xml('system-smirnoff.xml', reference_system)
            write_xml('openmm-smirnoff.xml', template_generated_system)
            raise Exception(f'Relative force deviation for {molecule.to_smiles()} ({relative_force_deviation}) exceeds threshold ({RELATIVE_FORCE_DEVIATION_TOLERANCE})')

class TestSMIRNOFFTemplateGenerator(TestGAFFTemplateGenerator):
    TEMPLATE_GENERATOR = SMIRNOFFTemplateGenerator

    def propagate_dynamics(self, molecule, system):
        """Run a few steps of dynamics to generate a perturbed configuration.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            molecule.conformers[0] is used as initial positions
        system : openmm.System
            System object for dynamics

        Returns
        -------
        new_molecule : openff.toolkit.topology.Molecule
            new_molecule.conformers[0] has updated positions

        """
        # Run some dynamics
        from openff.units.openmm import ensure_quantity
        from openmm import unit

        uses_old_api = hasattr(molecule.atoms[0], "element")

        temperature = 300 * unit.kelvin
        collision_rate = 1.0 / unit.picoseconds
        timestep = 1.0 * unit.femtoseconds
        nsteps = 100
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)
        context.setPositions(ensure_quantity(molecule.conformers[0], "openmm"))
        integrator.step(nsteps)
        # Copy the molecule, storing new conformer
        new_molecule = copy.deepcopy(molecule)
        new_positions = context.getState(getPositions=True).getPositions()

        new_molecule.conformers[0] = ensure_quantity(
            new_positions,
            "openmm" if uses_old_api else "openff",
        )


        del context, integrator

        return new_molecule

    def test_INSTALLED_FORCEFIELDS(self):
        """Test INSTALLED_FORCEFIELDS contains expected force fields"""
        expected_force_fields = [
            'openff-1.1.0',
            'openff-2.0.0',
            'smirnoff99Frosst-1.1.0',
        ]
        forbidden_force_fields = [
            'openff_unconstrained',
            'ff14sb_0.0.3',
        ]

        for expected in expected_force_fields:
            assert expected in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS

        for forbidden in forbidden_force_fields:
            assert forbidden not in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS

    def test_energies(self):
        """Test potential energies match between openff-toolkit and OpenMM ForceField"""

        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            if "ff14sb" in small_molecule_forcefield:
                continue
            if "tip" in small_molecule_forcefield:
                continue
            if "opc" in small_molecule_forcefield:
                continue

            # We cannot test openff-2.0.0-rc.1 because it triggers an openmm.OpenMMException due to an equilibrium angle > \pi
            # See https://github.com/openmm/openmm/issues/3185
            if "openff-2.0.0-rc.1" not in small_molecule_forcefield:
                continue

            print(f'Testing energies for {small_molecule_forcefield}...')
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=self.molecules, forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            for molecule in self.molecules:
                # Create OpenMM System using OpenMM app
                openmm_system = openmm_forcefield.createSystem(molecule.to_topology().to_openmm(), removeCMMotion=False, nonbondedMethod=NoCutoff)

                # Retrieve System generated by the SMIRNOFF typing engine
                smirnoff_system = generator.get_openmm_system(molecule)

                # Compare energies and forces
                self.compare_energies(molecule, openmm_system, smirnoff_system)

                # Run some dynamics
                molecule = self.propagate_dynamics(molecule, smirnoff_system)

                # Compare energies again
                self.compare_energies(molecule, openmm_system, smirnoff_system)


    def test_partial_charges_are_none(self):
        """Test parameterizing a small molecule with `partial_charges=None` instead
        of zeros (happens frequently in OFFTK>=0.7.0)"""
        molecule = Molecule.from_smiles('C=O')
        molecule.generate_conformers(n_conformers=1)
        #molecule._partial_charges = None
        assert (molecule.partial_charges is None) or np.all(molecule.partial_charges / unit.elementary_charge == 0)
        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            if "ff14sb" in small_molecule_forcefield:
                continue
            if "tip" in small_molecule_forcefield:
                continue
            if "opc" in small_molecule_forcefield:
                continue

            print(f'Testing energies for {small_molecule_forcefield}...')
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=[molecule], forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Create OpenMM System using OpenMM app
            openmm_system = openmm_forcefield.createSystem(molecule.to_topology().to_openmm(), removeCMMotion=False, nonbondedMethod=NoCutoff)
            smirnoff_system = generator.get_openmm_system(molecule)

    def test_version(self):
        """Test version"""
        # This test does not appear to test the version of anything in particular, but it fails sometimes
        # because old versions of the toolkit can't bring in new versions of some water models
        for forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator = SMIRNOFFTemplateGenerator(forcefield=forcefield)
            assert generator.forcefield == forcefield
            assert generator.smirnoff_filename.endswith(forcefield + '.offxml')
            assert os.path.exists(generator.smirnoff_filename)

    def test_bespoke_force_field(self):
        """
        Make sure a molecule can be parameterised using a bespoke force field passed as a string to
        the template generator.
        """

        custom_sage = OFFForceField("openff-2.0.0.offxml")
        # Create a simple molecule with one bond type
        ethane = Molecule.from_smiles("C")
        # Label ethane to get the bond type (not hard coded incase this changes in future)
        bond_parameter = custom_sage.label_molecules(ethane.to_topology())[0]["Bonds"][(0, 1)]
        # Edit the bond parameter
        bonds = custom_sage.get_parameter_handler("Bonds")
        new_parameter = bonds[bond_parameter.smirks]
        new_parameter.length = 2 * OFFUnit.angstrom

        # Use the custom sage passed as string to build a template and an openmm system
        generator = SMIRNOFFTemplateGenerator(molecules=ethane, forcefield=custom_sage.to_string())

        # Create a ForceField
        openmm_forcefield = openmm.app.ForceField()
        # Register the template generator
        openmm_forcefield.registerTemplateGenerator(generator.generator)
        # Use OpenMM app to generate the system
        openmm_system = openmm_forcefield.createSystem(ethane.to_topology().to_openmm(), removeCMMotion=False,
                                                       nonbondedMethod=NoCutoff)

        # Check the bond length has been updated
        forces = {force.__class__.__name__: force for force in openmm_system.getForces()}
        bond_force = forces["HarmonicBondForce"]
        for i in range(bond_force.getNumBonds()):
            _, _, length, _ = bond_force.getBondParameters(i)
            assert pytest.approx(length.value_in_unit(openmm.unit.angstrom)) == 2

    def test_14_scaling_from_offxml(self):
        """
        Make sure we can read lj14scale and coulomb14scale from the SMIRNOFF force field
        """
        custom_sage = OFFForceField("openff-2.0.0.offxml")
        custom_sage.get_parameter_handler("vdW").scale14 = 0.0
        custom_sage.get_parameter_handler("Electrostatics").scale14 = 0.0
        # Create a simplest 1-4 bond molecule
        ethane = Molecule.from_smiles("CC")

        # Use the custom sage passed as string to build a template and an openmm system
        generator = SMIRNOFFTemplateGenerator(molecules=ethane, forcefield=custom_sage.to_string())

        # Create a ForceField
        openmm_forcefield = openmm.app.ForceField()
        # Register the template generator
        openmm_forcefield.registerTemplateGenerator(generator.generator)
        # Use OpenMM app to generate the system
        openmm_system = openmm_forcefield.createSystem(
            ethane.to_topology().to_openmm(),
            removeCMMotion=False,
            nonbondedMethod=NoCutoff
        )

        # Find Nonbondedforce
        nb_force = [force for force in openmm_system.getForces() if force.__class__.__name__ == "NonbondedForce"][0]

        # Check all exceptions have q/eps == 0.
        exceptions = [nb_force.getExceptionParameters(i) for i in range(nb_force.getNumExceptions())]
        for exception in exceptions:
            assert exception[2]._value == 0.
            assert exception[4]._value == 0.

@pytest.mark.espaloma
class TestEspalomaTemplateGenerator(TestGAFFTemplateGenerator):
    TEMPLATE_GENERATOR = EspalomaTemplateGenerator

    def propagate_dynamics(self, molecule, system):
        """Run a few steps of dynamics to generate a perturbed configuration.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            molecule.conformers[0] is used as initial positions
        system : openmm.System
            System object for dynamics

        Returns
        -------
        new_molecule : openff.toolkit.topology.Molecule
            new_molecule.conformers[0] has updated positions

        """
        # Run some dynamics
        from openff.units.openmm import ensure_quantity
        from openmm import unit

        uses_old_api = hasattr(molecule.atoms[0], "element")

        temperature = 300 * unit.kelvin
        collision_rate = 1.0 / unit.picoseconds
        timestep = 1.0 * unit.femtoseconds
        nsteps = 100
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)
        context.setPositions(ensure_quantity(molecule.conformers[0], "openmm"))
        integrator.step(nsteps)
        # Copy the molecule, storing new conformer
        new_molecule = copy.deepcopy(molecule)
        new_positions = context.getState(getPositions=True).getPositions()

        new_molecule.conformers[0] = ensure_quantity(
            new_positions,
            "openmm" if uses_old_api else "openff",
        )

        # Clean up
        del context, integrator

        return new_molecule

    def test_retrieve_forcefields(self):
        """Test a force field can be retrieved"""
        # Test loading model by specifying version number
        generator = EspalomaTemplateGenerator(forcefield='espaloma-0.3.2')
        del generator
        # Test loading model from remote URL
        url = 'https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt'
        generator = EspalomaTemplateGenerator(forcefield=url)
        del generator
        # Test loading model from filename
        with tempfile.TemporaryDirectory() as tmpdirname:
            import urllib
            filename = os.path.join(tmpdirname, 'model.pt')
            urllib.request.urlretrieve(url, filename=filename)
            # Create a new database file
            generator = EspalomaTemplateGenerator(forcefield=filename)
            del generator

    def test_energies(self):
        """Test potential energies match between openff-toolkit and OpenMM ForceField"""

        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in EspalomaTemplateGenerator.INSTALLED_FORCEFIELDS:
            print(f'Testing energies for {small_molecule_forcefield}...')
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = EspalomaTemplateGenerator(molecules=self.molecules, forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            for molecule in self.molecules:
                # Create OpenMM System using OpenMM app
                openmm_system = openmm_forcefield.createSystem(molecule.to_topology().to_openmm(), removeCMMotion=False, nonbondedMethod=NoCutoff)

                # Retrieve System generated by Espaloma
                espaloma_system = generator.get_openmm_system(molecule)

                # Compare energies and forces
                self.compare_energies(molecule, openmm_system, espaloma_system)

                # Run some dynamics
                molecule = self.propagate_dynamics(molecule, espaloma_system)

                # Compare energies again
                self.compare_energies(molecule, openmm_system, espaloma_system)

    def test_partial_charges_are_none(self):
        """Test parameterizing a small molecule with `partial_charges=None` instead
        of zeros (happens frequently in OFFTK>=0.7.0)"""
        molecule = Molecule.from_smiles('C=O')
        molecule.generate_conformers(n_conformers=1)
        #molecule._partial_charges = None
        assert (molecule.partial_charges is None) or np.all(molecule.partial_charges / unit.elementary_charge == 0)
        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in EspalomaTemplateGenerator.INSTALLED_FORCEFIELDS:
            print(f'Testing energies for {small_molecule_forcefield}...')
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = EspalomaTemplateGenerator(molecules=[molecule], forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Create OpenMM System using OpenMM app
            openmm_system = openmm_forcefield.createSystem(molecule.to_topology().to_openmm(), removeCMMotion=False, nonbondedMethod=NoCutoff)
            smirnoff_system = generator.get_openmm_system(molecule)

    # def test_template_generator_kwargs(self):
    #     """Test """

    def test_keyword_arguments_default(self):
        """
        Test the default behavior for the keyword arguments, which is using "openff_unconstrained-2.0.0"
        as the reference forcefield and "from-molecule" as the charge method.

        Check charges are the same as the input, after passing through the template generator.
        """
        molecule = Molecule.from_smiles("C=O")
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges("am1bcc")  # Assign partial charges with off toolkit am1bcc method
        generator = EspalomaTemplateGenerator(molecules=[molecule], forcefield="espaloma-0.3.2")
        # Create forcefield object
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)
        # Create system
        system = forcefield.createSystem(molecule.to_topology().to_openmm(), nonbondedMethod=NoCutoff)
        # Make sure passing through the EspalomaGenerator didn't change the charges
        assert self.charges_are_equal(system, molecule), "Expected equal charges."
        # Assert the reference forcefield is the default "openff_unconstrained-2.0.0"
        default_ref_ff = "openff_unconstrained-2.0.0"
        generator_ref_ff = generator._reference_forcefield
        assert generator_ref_ff == default_ref_ff, f"Expected {default_ref_ff}, received {generator_ref_ff}."


    def test_keyword_arguments(self):
        """
        Test passing a custom keyword arguments dictionary to the generator behaves correctly.

        Specifically changing the reference forcefield to openff-2.1.0 and the charge method to nn.
        Which are not the default values.

        Check generated charges are different from the ones included in the input molecule.

        Returns
        -------

        """
        molecule = Molecule.from_smiles("C=O")
        molecule.generate_conformers(n_conformers=1)
        molecule.assign_partial_charges("am1bcc")  # Assign partial charges with off toolkit am1bcc method
        # Custom generator kwargs
        espaloma_generator_kwargs = {
            "reference_forcefield": "openff_unconstrained-2.1.0",
            "charge_method": "nn",
        }
        generator = EspalomaTemplateGenerator(molecules=[molecule],
                                              forcefield="espaloma-0.3.2",
                                              template_generator_kwargs=espaloma_generator_kwargs)
        # Create forcefield object
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)
        # Create system
        system = forcefield.createSystem(molecule.to_topology().to_openmm(), nonbondedMethod=NoCutoff)
        # Make sure passing through the EspalomaGenerator changes the charges
        assert not self.charges_are_equal(system, molecule), "Expected different charges"
        # Assert the reference forcefield is "openff_unconstrained-2.1.0"
        expected_ref_ff = "openff_unconstrained-2.1.0"
        generator_ref_ff = generator._reference_forcefield
        assert generator_ref_ff == expected_ref_ff, f"Expected {expected_ref_ff}, received {generator_ref_ff}."