import os, sys
import unittest
import tempfile

################################################################################
# Tests
################################################################################

class TestGAFFTemplateGenerator(unittest.TestCase):
    def setUp(self):
        # Read test molecules
        from openforcefield.topology import Molecule
        from openmmforcefields.utils import get_data_filename
        filename = get_data_filename("minidrugbank/MiniDrugBank.sdf")
        molecules = Molecule.from_file(filename, allow_undefined_stereo=True)
        # Select some small molecules for fast testing
        MAX_ATOMS = 24
        MAX_MOLECULES = 5
        molecules = [ molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS ]
        molecules = molecules[:MAX_MOLECULES]
        # Store molecules
        self.molecules = molecules

        # Suppress DEBUG logging from various packages
        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

    def test_create(self):
        """Test creation of an GAFFTemplateGenerator"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        # Create an empty generator
        generator = GAFFTemplateGenerator()
        # Create a generator that knows about a few molecules
        generator = GAFFTemplateGenerator(molecules=self.molecules)
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a new database file
            generator = GAFFTemplateGenerator(molecules=self.molecules, cache=cache)
            del generator
            # Reopen it (with cache still empty)
            generator = GAFFTemplateGenerator(molecules=self.molecules, cache=cache)
            del generator

    def test_parameterize(self):
        """Test parameterizing molecules with GAFFTemplateGenerator for all supported GAFF versions"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        # Test all supported GAFF versions
        for gaff_version in GAFFTemplateGenerator.SUPPORTED_GAFF_VERSIONS:
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = GAFFTemplateGenerator(molecules=self.molecules, gaff_version=gaff_version)
            # Create a ForceField with the appropriate GAFF version
            from simtk.openmm.app import ForceField
            forcefield = ForceField(generator.gaff_xml_filename)
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

    def test_add_molecules(self):
        """Test that molecules can be added to GAFFTemplateGenerator after its creation"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        # Create a generator that does not know about any molecules
        generator = GAFFTemplateGenerator()
        # Create a ForceField
        from simtk.openmm.app import ForceField
        forcefield = ForceField(generator.gaff_xml_filename)
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
        system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
        assert system.getNumParticles() == molecule.n_atoms

        # Add multiple molecules, including repeats
        generator.add_molecules(self.molecules)

        # Ensure all molecules can be parameterized
        print('5')
        for molecule in self.molecules:
            openmm_topology = molecule.to_topology().to_openmm()
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            assert system.getNumParticles() == molecule.n_atoms

    def test_cache(self):
        """Test cache capability of GAFFTemplateGenerator"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        from simtk.openmm.app import ForceField, NoCutoff
        gaff_version = '2.11'
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a generator that also has a database cache
            cache = os.path.join(tmpdirname, 'db.json')
            generator = GAFFTemplateGenerator(molecules=self.molecules, cache=cache, gaff_version=gaff_version)
            # Create a ForceField
            forcefield = ForceField(generator.gaff_xml_filename)
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
                generator : GAFFTemplateGenerator
                    The generator whose cache should be examined
                n_expected : int
                    Number of expected records
                """
                from tinydb import TinyDB
                db = TinyDB(generator._cache)
                table = db.table(generator.gaff_version)
                db_entries = table.all()
                db.close()
                n_entries = len(db_entries)
                assert (n_entries, n_expected), \
                    "Expected {} entries but database has {}\n db contents: {}".format(n_expected, n_entries, db_entries)

            check_cache(generator, len(self.molecules))

            # Clean up, forcing closure of database
            del forcefield, generator

            # Create a generator that also uses the database cache but has no molecules
            print('Creating new generator with just cache...')
            generator = GAFFTemplateGenerator(cache=cache, gaff_version=gaff_version)
            # Check database still contains the molecules we expect
            check_cache(generator, len(self.molecules))
            # Create a ForceField
            forcefield = ForceField(generator.gaff_xml_filename)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules; this should succeed
            for molecule in self.molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)

            # Changing the GAFF version should use a different table with no molecule entries
            generator = GAFFTemplateGenerator(cache=cache, gaff_version='1.81')
            # Check database contains no molecules
            check_cache(generator, 0)
