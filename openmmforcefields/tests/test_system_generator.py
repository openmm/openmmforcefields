import os, sys
import unittest
import tempfile

from openmmforcefields.utils import get_data_filename, Timer

################################################################################
# Tests
################################################################################

class TestSystemGenerator(unittest.TestCase):
    def setUp(self):
        # Create test topologies
        self.topologies = list()

        # Create protein-ligand system
        from openforcefield.topology import Molecule
        ligand_off_molecule = Molecule(get_data_filename("BRD4/sdf/ligand-1.sdf"))

        # Read in the coordinates of the ligand from the PDB file
        from simtk.openmm.app import PDBFile
        ligand_pdbfile = PDBFile('ligand.pdb')

        # Com

        # TODO: Create other test topologies
        # TODO: Protein-only
        # TODO: Protein-ligand topology
        # TODO: Solvated protein-ligand topology
        # TODO: Host-guest topology

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

class TestGAFFSystemGenerator(TestSystemGenerator):

    def test_create(self):
        """Test GAFFSystemGenerator creation"""
        from openmmforcefields.generators import GAFFSystemGenerator
        # Create an empty system generator
        generator = GAFFSystemGenerator()
        # Create a generator that defines a protein force field and GAFF version
        generator = GAFFSystemGenerator(forcefields=['amber/protein.ff14SB', 'amber/tip3p.xml'], gaff_version='2.11')
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a new database file
            generator = GAFFSystemGenerator(cache=cache)
            del generator
            # Reopen it (with cache still empty)
            generator = GAFFSystemGenerator(cache=cache)
            del generator

    def test_parameterize(self):
        """Test parameterizing molecules with GAFFTemplateGenerator"""
        from openmmforcefields.generators import GAFFTemplateGenerator, GAFFSystemGenerator
        for gaff_version in GAFFTemplateGenerator.SUPPORTED_GAFF_VERSIONS:
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            forcefields = ['amber/protein.ff14SB', 'tip3p.xml']
            generator = GAFFSystemGenerator(forcefields=forcefields, gaff_version=gaff_version)
            # Parameterize some systems
            for topology in self.topologies:
                with Timer() as t1:
                    system = generator.create_system(topology)
                assert system.getNumParticles() == topology.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = generator.create_system(topology)
                assert system.getNumParticles() == topology.n_atoms
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
