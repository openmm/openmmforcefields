import os, sys
import unittest
import tempfile

from openmmforcefields.utils import get_data_filename, Timer

################################################################################
# Tests
################################################################################

class TestGAFFTemplateGenerator(unittest.TestCase):
    def setUp(self):
        # Read test molecules
        from openforcefield.topology import Molecule
        filename = get_data_filename("minidrugbank/MiniDrugBank.sdf")
        molecules = Molecule.from_file(filename, allow_undefined_stereo=True)
        # Only use first file molecules for testing
        molecules = molecules[:5]
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
            for molecule in self.molecules:
                topology = molecule.to_openmm_topology()
                with Timer() as t1:
                    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == molecule.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = forcefield.createSystem(topology, nonbondedMethod=NoCutoff)
                assert system.getNumParticles() == molecule.n_atoms
                assert (t2.interval < t1.interval)

    def test_add_molecules(self):
        """Test that OEMols can be added to GAFFTemplateGenerator later"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        # Create a generator that does not know about any molecules
        generator = GAFFTemplateGenerator(gaff_version='2.11')
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
            system = forcefield.createSystem(molecule.to_openmm_topology(), nonbondedMethod=NoCutoff)
        except ValueError as e:
            # Exception 'No template found...' is expected
            assert str(e).startswith('No template found')

        # Now add the molecule to the generator and ensure parameterization passes
        generator.add_molecules(molecule)
        system = forcefield.createSystem(molecule.to_openmm_topology(), nonbondedMethod=NoCutoff)
        assert system.getNumParticles() == molecule.n_atoms

        # Add multiple molecules, including repeats
        generator.add_molecules(self.molecules)

        # Ensure all molecules can be parameterized
        for molecule in self.molecules:
            system = forcefield.createSystem(molecule.to_openmm_topology(), nonbondedMethod=NoCutoff)
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
            for oemol in self.oemols:
                forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)

            # Check database contents
            def check_database(generator):
                from tinydb import TinyDB
                db = TinyDB(generator._cache)
                db_entries = db.all()
                db.close()
                nentries = len(db_entries)
                nmolecules = len(self.oemols)
                assert (nmolecules == nentries), \
                    "Expected {} entries but database only has {}\n db contents: {}".format(nmolecules, nentries, db_entries)

            check_database(generator)
            # Clean up, forcing closure of database
            del forcefield, generator

            # Create a generator that also uses the database cache but has no molecules
            print('Creating new generator with just cache...')
            generator = GAFFTemplateGenerator(cache=cache, gaff_version=gaff_version)
            # Check database contents
            check_database(generator)
            # Create a ForceField
            forcefield = ForceField(generator.gaff_xml_filename)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize the molecules; this should succeed
            for oemol in self.oemols:
                from openeye import oechem
                forcefield.createSystem(generateTopologyFromOEMol(oemol), nonbondedMethod=NoCutoff)

            # Check that using an incompatible cache will fail
            from openmmforcefields.generators import IncompatibleGAFFVersion
            try:
                generator = GAFFTemplateGenerator(cache=cache, gaff_version='1.81')
            except IncompatibleGAFFVersion as e:
                pass
