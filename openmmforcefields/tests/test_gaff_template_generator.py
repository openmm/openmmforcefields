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
                assert (n_entries == n_expected), \
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

    def test_add_solvent(self):
        """Test using simtk.opnmm.app.Modeller to add solvent to a small molecule parameterized by GAFFTemplateGenerator"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        gaff_version = '2.11'
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
        generator = GAFFTemplateGenerator(molecules=self.molecules, gaff_version=gaff_version)
        # Add to the forcefield object
        forcefield.loadFile(generator.gaff_xml_filename)
        forcefield.registerTemplateGenerator(generator.generator)
        # Add solvent to a system containing a small molecule
        # This should succeed
        modeller.addSolvent(forcefield, model='tip3p', padding=6.0*unit.angstroms)

    def test_jacs_ligands(self):
        """Use GAFFTemplateGenerator to parameterize the Schrodinger JACS set of ligands"""
        from simtk.openmm.app import ForceField, NoCutoff
        gaff_version = '2.11'

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
            print(f'Read {len(molecules)} molecules from {sdf_filename}')

            # Create GAFF template generator with local cache
            cache_filename = os.path.join(get_data_filename(os.path.join('perses_jacs_systems', system_name)), 'cache.json')
            from openmmforcefields.generators import GAFFTemplateGenerator
            generator = GAFFTemplateGenerator(molecules=molecules, gaff_version=gaff_version, cache=cache_filename)

            # Create a ForceField
            forcefield = ForceField(generator.gaff_xml_filename)
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
