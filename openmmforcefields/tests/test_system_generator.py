import os, sys
import unittest
import tempfile

from openmmforcefields.utils import get_data_filename, Timer

################################################################################
# Tests
################################################################################

class TestSystemGenerator(unittest.TestCase):
    def setUp(self):
        self.testsystems = dict()
        for (system_name, prefix) in [
            ('bace', 'Bace'),
            ('cdk1', 'CDK2'),
            ('jnk1', 'Jnk1'),
            ('mcl1', 'MCL1'),
            ('p38', 'p38'),
            ('ptp1b', 'PTP1B'),
            ('thrombin', 'Thrombin'),
            ('tyk2', 'Tyk2'),
        ]:
            # Load protein
            from simtk.openmm.app import PDBFile
            pdb_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_protein_fixed.pdb'))
            pdbfile = PDBFile(pdb_filename)

            # Load molecules
            from openforcefield.topology import Molecule
            from openmmforcefields.utils import get_data_filename
            sdf_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_ligands.sdf'))
            molecules = Molecule.from_file(sdf_filename, allow_undefined_stereo=True)
            n_molecules = len(molecules)
            print(f'Read {n_molecules} molecules from {sdf_filename}')

            # Create ParmEd Structure objects
            print('Creating protein:ligand complexes')
            import parmed
            protein_structure = parmed.load_file(pdb_filename)
            molecules_structure = parmed.load_file(sdf_filename)
            complex_structures = [ (protein_structure + molecules_structure[index]) for index in range(n_molecules) ]

            # Store
            testsystem = {
                'name' : system_name,
                'protein_pdbfile' : pdbfile,
                'molecules' : molecules,
                'complex_structures' : complex_structures
                }
            self.testsystems[system_name] = testsystem

        # TODO: Create other test topologies
        # TODO: Protein-only
        # TODO: Protein-ligand topology
        # TODO: Solvated protein-ligand topology
        # TODO: Host-guest topology

        # Select AMBER force fields
        self.amber_forcefields = ['amber/protein.ff14SB', 'amber/tip3p_standard.xml']

        # Suppress DEBUG logging from various packages
        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

    def test_create(self):
        """Test GAFFSystemGenerator creation"""
        # Create an empty system generator
        generator = self.SystemGenerator()
        # Create a generator that defines a protein force field and GAFF version
        generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                        **self.extra_generator_arguments)
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, 'db.json')
            # Create a new database file
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                            cache=cache, **self.extra_generator_arguments)
            del generator
            # Reopen it (with cache still empty)
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                            cache=cache, **self.extra_generator_arguments)
            del generator

    def test_parameterize_molecules(self):
        """Test parameterizing molecules in vacuum with GAFFTemplateGenerator"""
        from openmmforcefields.generators import GAFFSystemGenerator
        for gaff_version in GAFFTemplateGenerator.SUPPORTED_GAFF_VERSIONS:
            # Create a generator that knows about a few molecules
            testsystem = self.testsystems['mcl1']
            molecules = testsystem['molecules']
            from simtk.openmm.app import NoCutoff
            forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                             forcefield_kwargs=forcefield_kwargs,
                                             molecules=molecules, **self.extra_generator_arguments)
            # Parameterize some molecules
            from openmmforcefields.utils import Timer
            for molecule in molecules[:3]:
                openmm_topology = molecule.to_topology().to_openmm()
                with Timer() as t1:
                    system = generator.create_system(openmm_topology)
                assert system.getNumParticles() == molecule.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = generator.create_system(topology)
                assert system.getNumParticles() == molecule.n_atoms
                assert (t2.interval < t1.interval)

    def test_add_molecules(self):
        """Test that Molecules can be added to GAFFSystemGenerator later"""
        from openmmforcefields.generators import GAFFTemplateGenerator
        # Create a generator that does not know about any molecules
        testsystem = self.testsystems['mcl1']
        molecules = testsystem['molecules']
        from simtk.openmm.app import NoCutoff
        forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }
        generator = self.SystemGenerator(forcefields=self.amber_forcefields,
            forcefield_kwargs=forcefield_kwargs, **self.extra_generator_arguments)

        # Check that parameterizing a molecule fails
        molecule = molecules[0]
        openmm_topology = molecule.to_topology().to_openmm()
        try:
            # This should fail with an exception
            system = forcefield.create_system(openmm_topology)
        except ValueError as e:
            # Exception 'No template found...' is expected
            assert str(e).startswith('No template found')

        # Now add the molecule to the generator and ensure parameterization passes
        generator.add_molecules(molecule)
        system = forcefield.createSystem(openmm_topology)
        assert system.getNumParticles() == molecule.n_atoms

        # Add multiple molecules, including repeats
        generator.add_molecules(molecules)

        # Ensure all molecules can be parameterized
        for molecule in molecules[:3]:
            system = forcefield.createSystem(openmm_topology)
            assert system.getNumParticles() == molecule.n_atoms

    def test_complex(self):
        """Test parameterizing a protein:ligand complex in vacuum"""
        for testsystem in self.testsystems:
            print(f'Testing parameterization of {testsystem.name}')
            molecules = testsystem['molecules']
            from simtk.openmm.app import NoCutoff
            forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                forcefield_kwargs=forcefield_kwargs,
                molecules=molecules, **self.extra_generator_arguments)
            # Parameterize a complex from the set
            complex_structure = testsystem['complex_structures'][0]
            openmm_topology = complex_structure.topology
            system = generator.create_system(openmm_topology)
            assert system.getNumParticles() == len(complex_structure.atoms)

    # TODO: Test parameterization of a protein:ligand complex in solvent

class TestGAFFSystemGenerator(TestSystemGenerator):
    """Test the GAFFSystemGenerator"""
    from openmmforcefields.generators import GAFFSystemGenerator
    SystemGenerator = GAFFSystemGenerator
    extra_generator_arguments = { 'gaff_version' : '2.11' }

class TestSMIRNOFFSystemGenerator(TestSystemGenerator):
    """Test the SMIRNOFFSystemGenerator"""
    from openmmforcefields.generators import GAFFSystemGenerator
    SystemGenerator = GAFFSystemGenerator
    extra_generator_arguments = { 'smirnoff_filename' : 'openforcefield-1.0.0' }
