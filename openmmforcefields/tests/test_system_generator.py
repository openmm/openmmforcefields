import os, sys
import unittest
import tempfile

from openmmforcefields.utils import get_data_filename, Timer

################################################################################
# Tests
################################################################################

class TestSystemGenerator(object):

    SMALL_MOLECULE_FORCEFIELDS = ['gaff-1.81', 'gaff-2.11', 'smirnoff99Frosst-1.1.0', 'openforcefield-1.0.0']

    # AMBER force field combination to test
    amber_forcefields = ['amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml']

    def filter_molecules(self, molecules):
        """
        Filter molecules to speed up tests, especially on travis.

        Parameters
        ----------
        molecules : list of openforcefield.topology.Molecule
            The input list of molecules to be filtered

        Returns
        -------
        molecules : list of openforcefield.topology.Molecule
            The filtered list of molecules to be filtered

        """
        # TODO: Eliminate molecules without fully-specified stereochemistry
        # Select some small molecules for fast testing
        MAX_ATOMS = 24
        molecules = [ molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS ]
        # Cut down number of tests for travis
        import os
        MAX_MOLECULES = 20
        if 'TRAVIS' in os.environ:
            MAX_MOLECULES = 3
        molecules = molecules[:MAX_MOLECULES]

        return molecules

    """Base class for SystemGenerator tests."""
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
            print(f'Read {len(molecules)} molecules from {sdf_filename}')
            # Filter molecules as appropriate
            molecules = self.filter_molecules(molecules)
            print(f'{len(molecules)} molecules remain after filtering')

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
        # Suppress DEBUG logging from various packages

        import logging
        for name in ['parmed', 'matplotlib']:
            logging.getLogger(name).setLevel(logging.WARNING)

    def test_create(self):
        """Test SystemGenerator creation with only OpenMM ffxml force fields"""
        from openmmforcefields.generators import SystemGenerator

        # Create an empty system generator
        generator = self.SystemGenerator()

    def test_barostat(self):
        """Test that barostat addition works correctly"""
        # Create a protein SystemGenerator
        from openmmforcefields.generators import SystemGenerator
        generator = self.SystemGenerator(forcefields=self.amber_forcefields)

        # Create a template barostat
        from simtk.openmm import MonteCarloBarostat
        from simtk import unit
        pressure = 0.95 * unit.atmospheres
        temperature = 301.0 * unit.kelvin
        frequency = 23
        generator.barostat = MonteCarloBarostat(pressure, temperature, frequency)

        # Load a PDB file
        import os
        from simtk.openmm.app import PDBFile
        pdb_filename = get_data_filename(os.path.join('perses_jacs_systems', 'mcl1', 'MCL1_protein_fixed.pdb'))
        pdbfile = PDBFile(pdb_filename)

        # Create a System
        system = generator.create_system(pdbfile.topology)

        # Check barostat is present
        forces = { force.__class__.__name__ : force for force in system.getForces() }
        assert 'MonteCarloBarostat' in forces.keys()

        # Check barostat parameters
        force = forces['MonteCarloBarostat']
        assert force.getDefaultPressure() == pressure
        assert force.getDefaultTemperature() == temperature
        assert force.getFrequency() == frequency

    def test_create_with_template_generator(self):
        """Test SystemGenerator creation with small molecule residue template generators"""
        from openmmforcefields.generators import SystemGenerator

        for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
            # Create a generator that defines AMBER and small molecule force fields
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield)

            # Create a generator that also has a database cache
            with tempfile.TemporaryDirectory() as tmpdirname:
                cache = os.path.join(tmpdirname, 'db.json')
                # Create a new database file
                generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                    cache=cache, small_molecule_forcefield=small_molecule_forcefield)
                del generator
                # Reopen it (with cache still empty)
                generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                    cache=cache, small_molecule_forcefield=small_molecule_forcefield)
                del generator

    def test_parameterize_molecules_from_creation(self):
        """Test that SystemGenerator can parameterize pre-specified molecules in vacuum"""
        from openmmforcefields.generators import SystemGenerator
        from simtk.openmm.app import NoCutoff
        forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }

        for testsystem in self.testsystems:
            molecules = testsystem['molecules']

            for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator for this force field
                generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
                                                forcefield_kwargs=forcefield_kwargs,
                                                molecules=molecules)

                # Parameterize molecules
                from openmmforcefields.utils import Timer
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    with Timer() as t1:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    # Molecule should now be cached
                    with Timer() as t2:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    assert (t2.interval < t1.interval)

    def test_parameterize_molecules_specified_during_create_system(self):
        """Test that SystemGenerator can parameterize molecules specified during create_system"""
        from openmmforcefields.generators import SystemGenerator
        from simtk.openmm.app import NoCutoff
        forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }

        for testsystem in self.testsystems:
            molecules = testsystem['molecules']

            for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator for this force field
                generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
                                                forcefield_kwargs=forcefield_kwargs)

                # Parameterize molecules
                from openmmforcefields.utils import Timer
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    # Specify molecules during system creation
                    system = generator.create_system(openmm_topology, molecules=molecules)

    def test_add_molecules(self):
        """Test that Molecules can be added to SystemGenerator later"""
        from openmmforcefields.generators import SystemGenerator
        from simtk.openmm.app import NoCutoff
        forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }

        for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
            # Create a SystemGenerator for this force field
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                            small_molecule_forcefield=small_molecule_forcefield,
                                            forcefield_kwargs=forcefield_kwargs,
                                            molecules=molecules)

            # Add molecules for each test system separately
            for testsystem in self.testsystems:
                molecules = testsystem['molecules']
                # Add molecules
                generator.add_molecules(molecules)

                # Parameterize molecules
                from openmmforcefields.utils import Timer
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    with Timer() as t1:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    # Molecule should now be cached
                    with Timer() as t2:
                        system = generator.create_system(topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    assert (t2.interval < t1.interval)

def test_cache(self):
    """Test that SystemGenerator correctly manages a cache"""
    from openmmforcefields.generators import SystemGenerator
    from simtk.openmm.app import NoCutoff
    from openmmforcefields.utils import Timer
    forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }

    timing = dict() # timing[(small_molecule_forcefield, smiles)] is the time (in seconds) to parameterize molecule the first time
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Create a single shared cache for all force fields
        cache = os.path.join(tmpdirname, 'db.json')
        # Test that we can parameterize all molecules for all test systems
        for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
            # Create a SystemGenerator
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                            small_molecule_forcefield=small_molecule_forcefield,
                                            forcefield_kwargs=forcefield_kwargs,
                                            molecules=molecules,
                                            cache=cache)
            # Add molecules for each test system separately
            for testsystem in self.testsystems:
                molecules = testsystem['molecules']
                # Add molecules
                generator.add_molecules(molecules)

                # Parameterize molecules
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    with Timer() as elapsed_time:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    # Record time
                    timing[(small_molecule_forcefield, molecule.to_smiles())] = elapsed_time

        # Molecules should now be cached; test timing is faster the second time
        # Test that we can parameterize all molecules for all test systems
        for small_molecule_forcefield in SMALL_MOLECULE_FORCEFIELDS:
            # Create a SystemGenerator
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                                            small_molecule_forcefield=small_molecule_forcefield,
                                            forcefield_kwargs=forcefield_kwargs,
                                            molecules=molecules,
                                            cache=cache)
            # Add molecules for each test system separately
            for testsystem in self.testsystems:
                molecules = testsystem['molecules']
                # We don't need to add molecules that are already defined in the cache

                # Parameterize molecules
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    with Timer() as elapsed_time:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    # Check that this was faster than the first time
                    assert (elapsed_time < timing[(small_molecule_forcefield, molecule.to_smiles())])

    def test_complex(self):
        """Test parameterizing a protein:ligand complex in vacuum"""
        from openmmforcefields.generators import SystemGenerator
        for testsystem in self.testsystems:
            print(f'Testing parameterization of {testsystem.name} in vacuum')
            molecules = testsystem['molecules']
            from simtk.openmm.app import NoCutoff
            forcefield_kwargs = { 'nonbondedMethod' : NoCutoff }
            generator = self.SystemGenerator(forcefields=self.amber_forcefields,
                forcefield_kwargs=forcefield_kwargs,
                molecules=molecules)
            # Parameterize a complex from the set
            complex_structure = testsystem['complex_structures'][0]
            openmm_topology = complex_structure.topology
            system = generator.create_system(openmm_topology)
            assert system.getNumParticles() == len(complex_structure.atoms)

    # TODO: Test parameterization of a protein:ligand complex in solvent with ions
