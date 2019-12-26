import os, sys
import unittest
import tempfile

from openmmforcefields.generators import GAFFTemplateGenerator
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

################################################################################
# Tests
################################################################################

class TestGAFFTemplateGenerator(unittest.TestCase):
    TEMPLATE_GENERATOR = GAFFTemplateGenerator

    def setUp(self):
        # Read test molecules
        from openforcefield.topology import Molecule
        from openmmforcefields.utils import get_data_filename
        filename = get_data_filename("minidrugbank/MiniDrugBank-without-unspecified-stereochemistry.sdf")
        molecules = Molecule.from_file(filename, allow_undefined_stereo=True)
        # Select some small molecules for fast testing
        MAX_ATOMS = 20
        molecules = [ molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS ]
        # Cut down number of tests for travis
        import os
        if 'TRAVIS' in os.environ:
            MAX_MOLECULES = 2
            molecules = molecules[:MAX_MOLECULES]
        # Store molecules
        self.molecules = molecules

        # Suppress DEBUG logging from various packages
        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

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
        from simtk.openmm.app import ForceField
        forcefield = ForceField()
        # Register the template generator
        forcefield.registerTemplateGenerator(generator.generator)

        # Check that parameterizing a molecule fails
        molecule = self.molecules[0]
        from simtk.openmm.app import NoCutoff
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
            from simtk.openmm.app import PDBFile
            PDBFile.writeFile(openmm_topology, molecule.conformers[0])
            raise e
        assert system.getNumParticles() == molecule.n_atoms

        # Add multiple molecules, including repeats
        generator.add_molecules(self.molecules)

        # Ensure all molecules can be parameterized
        for molecule in self.molecules:
            openmm_topology = molecule.to_topology().to_openmm()
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert system.getNumParticles() == molecule.n_atoms

    def test_cache(self):
        """Test template generator cache capability"""
        from simtk.openmm.app import ForceField, NoCutoff
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
                    "Expected {} entries but database has {}\n db contents: {}".format(n_expected, n_entries, db_entries)

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
        """Test using simtk.opnmm.app.Modeller to add solvent to a small molecule parameterized by template generator"""
        # Select a molecule to add solvent around
        from simtk.openmm.app import NoCutoff, Modeller
        from simtk import unit
        molecule = self.molecules[0]
        openmm_topology = molecule.to_topology().to_openmm()
        openmm_positions = molecule.conformers[0]
        # Try adding solvent without residue template generator; this will fail
        from simtk.openmm.app import ForceField
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
        from simtk.openmm.app import ForceField, NoCutoff
        jacs_systems = {
            'bace'     : { 'ligand_sdf_filename' : 'Bace_ligands.sdf' },
            'cdk2'     : { 'ligand_sdf_filename' : 'CDK2_ligands.sdf' },
            'jnk1'     : { 'ligand_sdf_filename' : 'Jnk1_ligands.sdf' },
            'mcl1'     : { 'ligand_sdf_filename' : 'MCL1_ligands.sdf' },
            'p38'      : { 'ligand_sdf_filename' : 'p38_ligands.sdf' },
            'ptp1b'    : { 'ligand_sdf_filename' : 'PTP1B_ligands.sdf' },
            'thrombin' : { 'ligand_sdf_filename' : 'Thrombin_ligands.sdf' },
            'tyk2'     : { 'ligand_sdf_filename' : 'Tyk2_protein.pdb' },
        }
        for system_name in jacs_systems:
            # Load molecules
            ligand_sdf_filename = jacs_systems[system_name]['ligand_sdf_filename']
            print(f'Reading molecules from {ligand_sdf_filename} ...')
            from openforcefield.topology import Molecule
            from openmmforcefields.utils import get_data_filename
            sdf_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, ligand_sdf_filename))
            molecules = Molecule.from_file(sdf_filename, allow_undefined_stereo=True)
            # Ensure this is a list
            try:
                nmolecules = len(molecules)
            except TypeError:
                molecules = [molecules]

            print(f'Read {len(molecules)} molecules from {sdf_filename}')

            # Create GAFF template generator with local cache
            cache_filename = os.path.join(get_data_filename(os.path.join('perses_jacs_systems', system_name)), 'cache.json')
            generator = self.TEMPLATE_GENERATOR(molecules=molecules, cache=cache_filename)

            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)

            # Parameterize all molecules
            print(f'Caching all molecules for {system_name} at {cache_filename} ...')
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

    # TODO: Test JACS protein-ligand systems

    def test_parameterize(self):
        """Test parameterizing molecules with GAFFTemplateGenerator for all supported GAFF versions"""
        # Test all supported GAFF versions
        for gaff_version in GAFFTemplateGenerator.SUPPORTED_GAFF_VERSIONS:
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = GAFFTemplateGenerator(molecules=self.molecules, gaff_version=gaff_version)
            # Create a ForceField with the appropriate GAFF version
            from simtk.openmm.app import ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            from simtk.openmm.app import NoCutoff
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

class TestSMIRNOFFTemplateGenerator(TestGAFFTemplateGenerator):
    TEMPLATE_GENERATOR = SMIRNOFFTemplateGenerator

    def test_parameterize(self):
        """Test parameterizing molecules with SMIRNOFFTemplateGenerator for all installed SMIRNOFF force fields"""
        # Test all supported SMIRNOFF force fields
        for smirnoff in SMIRNOFFTemplateGenerator.INSTALLED_SMIRNOFF_FORCEFIELDS:
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=self.molecules, smirnoff=smirnoff)
            # Create a ForceField
            from simtk.openmm.app import ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            from simtk.openmm.app import NoCutoff
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

    @staticmethod
    def compute_energy(system, positions):
        """Compute potential energy and Force components for an OpenMM system.

        Parameters
        ----------
        system : simtk.openmm.System
            The System object
        positions : simtk.unit.Quantity of shape (nparticles,3) with units compatible with nanometers
            The positions for which energy is to be computed

        Returns
        -------
        openmm_energy : dict of str : simtk.unit.Quantity
            openmm_energy['total'] is the total potential energy
            openmm_energy['components'][forcename] is the potential energy for the specified component force
        """
        import copy
        system = copy.deepcopy(system)
        for index, force in enumerate(system.getForces()):
            force.setForceGroup(index)
        from simtk import openmm
        platform = openmm.Platform.getPlatformByName('Reference')
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        openmm_energy = {
            'total' : context.getState(getEnergy=True).getPotentialEnergy(),
            'components' : { system.getForce(index).__class__.__name__ : context.getState(getEnergy=True, groups=(1 << index)).getPotentialEnergy() for index in range(system.getNumForces()) },
            }
        del context, integrator
        return openmm_energy

    @classmethod
    def compare_energies(cls, molecule, off_forcefield, openmm_forcefield):
        """Compare energies between Open Force Field Initiative and OpenMM ForceField objects"""
        from simtk import openmm

        # Create OpenMM System using openfocefield
        smirnoff_system = off_forcefield.create_openmm_system(molecule.to_topology())
        smirnoff_energy = cls.compute_energy(smirnoff_system, molecule.conformers[0])
        with open('system-smirnoff.xml', 'w') as outfile:
            outfile.write(openmm.XmlSerializer.serialize(smirnoff_system))

        # Create OpenMM System using OpenMM app
        from simtk.openmm import app
        openmm_system = openmm_forcefield.createSystem(molecule.to_topology().to_openmm(), removeCMMotion=False, onbondedMethod=app.NoCutoff)
        openmm_energy = cls.compute_energy(openmm_system, molecule.conformers[0])
        with open('system-openmm.xml', 'w') as outfile:
            outfile.write(openmm.XmlSerializer.serialize(openmm_system))

        # TODO: Compare energies
        from simtk import unit
        ENERGY_DEVIATION_TOLERANCE = 1.0e-2 * unit.kilocalories_per_mole
        delta = (openmm_energy['total'] - smirnoff_energy['total'])
        if abs(delta) > ENERGY_DEVIATION_TOLERANCE:
            # Show breakdown by components
            print('Energy components:')
            print(f"{'component':24} {'OpenMM (kcal/mol)':20} {'SMIRNOFF (kcal/mol)':20}")
            for key in openmm_energy['components'].keys():
                openmm_component_energy = openmm_energy['components'][key]
                smirnoff_component_energy = smirnoff_energy['components'][key]
                print(f'{key:24} {(openmm_component_energy/unit.kilocalories_per_mole):20.3f} {(smirnoff_component_energy/unit.kilocalories_per_mole):20.3f} kcal/mol')
            print(f'{"TOTAL":24} {(openmm_energy["total"]/unit.kilocalories_per_mole):20.3f} {(smirnoff_energy["total"]/unit.kilocalories_per_mole):20.3f} kcal/mol')
            raise Exception(f'Energy deviation ({delta/unit.kilocalories_per_mole} kcal/mol) exceeds threshold ({ENERGY_DEVIATION_TOLERANCE})')

    def test_energies(self):
        """Test potential energies match between openforcefield and OpenMM ForceField"""
        # Test all supported SMIRNOFF force fields
        for smirnoff in SMIRNOFFTemplateGenerator.INSTALLED_SMIRNOFF_FORCEFIELDS:
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=self.molecules, smirnoff=smirnoff)
            # Create a ForceField
            import simtk
            openmm_forcefield = simtk.openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Create openforcefield ForceField object
            import openforcefield
            off_forcefield = openforcefield.typing.engines.smirnoff.ForceField(smirnoff)
            # Parameterize some molecules
            from simtk.openmm.app import NoCutoff
            for molecule in self.molecules:
                self.compare_energies(molecule, off_forcefield, openmm_forcefield)
