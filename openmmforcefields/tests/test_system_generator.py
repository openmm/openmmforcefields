import os, sys
import unittest
import tempfile

from openmmforcefields.utils import get_data_filename

from openmmforcefields.generators import SystemGenerator

################################################################################
# Tests
################################################################################

class TestSystemGenerator(unittest.TestCase):
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
        MAX_ATOMS = 40
        molecules = [ molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS ]
        # Cut down number of tests for travis
        import os
        MAX_MOLECULES = 10
        if 'TRAVIS' in os.environ:
            MAX_MOLECULES = 3
        molecules = molecules[:MAX_MOLECULES]

        return molecules

    """Base class for SystemGenerator tests."""
    def setUp(self):
        self.testsystems = dict()
        for (system_name, prefix) in [
            # TODO: Uncomment these after we fix input files
            #('bace', 'Bace'),
            #('cdk1', 'CDK2'),
            ('jnk1', 'Jnk1'),
            #('mcl1', 'MCL1'),
            #('p38', 'p38'),
            #('ptp1b', 'PTP1B'),
            #('thrombin', 'Thrombin'),
            #('tyk2', 'Tyk2'),
        ]:
            # Load protein
            from simtk.openmm.app import PDBFile
            pdb_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_protein_fixed.pdb'))
            pdbfile = PDBFile(pdb_filename)

            # Load molecules
            from openforcefield.topology import Molecule
            sdf_filename = get_data_filename(os.path.join('perses_jacs_systems', system_name, prefix + '_ligands.sdf'))

            molecules = Molecule.from_file(sdf_filename, allow_undefined_stereo=True)
            print(f'Read {len(molecules)} molecules from {sdf_filename}')
            # Filter molecules as appropriate
            molecules = self.filter_molecules(molecules)
            n_molecules = len(molecules)
            print(f'{n_molecules} molecules remain after filtering')
            if n_molecules == 0:
                continue

            # Create structures
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
        # Create an empty system generator
        generator = SystemGenerator()

    def test_barostat(self):
        """Test that barostat addition works correctly"""
        # Create a protein SystemGenerator
        generator = SystemGenerator(forcefields=self.amber_forcefields)

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
        pdb_filename = get_data_filename(os.path.join('perses_jacs_systems', 'bace', 'Bace_protein_fixed.pdb'))
        pdbfile = PDBFile(pdb_filename)

        # Delete hydrogens from terminal protein residues
        # TODO: Fix the input files so we don't need to do this
        from simtk.openmm import app
        modeller = app.Modeller(pdbfile.topology, pdbfile.positions)
        residues = [residue for residue in modeller.topology.residues() if residue.name != 'UNL']
        termini_ids = [residues[0].id, residues[-1].id]
        #hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol in ['H'] and atom.residue.name != 'UNL']
        hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol in ['H'] and atom.residue.id in termini_ids]
        modeller.delete(hs)
        from simtk.openmm.app import PDBFile
        modeller.addHydrogens()

        # Create a System
        system = generator.create_system(modeller.topology)

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
        for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
            # Create a generator that defines AMBER and small molecule force fields
            generator = SystemGenerator(forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield)

            # Create a generator that also has a database cache
            with tempfile.TemporaryDirectory() as tmpdirname:
                cache = os.path.join(tmpdirname, 'db.json')
                # Create a new database file
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                    cache=cache, small_molecule_forcefield=small_molecule_forcefield)
                del generator
                # Reopen it (with cache still empty)
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                    cache=cache, small_molecule_forcefield=small_molecule_forcefield)
                del generator

    def test_forcefield_kwargs(self):
        """Test that forcefield_kwargs and nonbonded method specifications work correctly"""
        from simtk import unit
        forcefield_kwargs = { 'hydrogenMass' : 4*unit.amu }
        from openmmforcefields.generators import SystemGenerator

        for name, testsystem in self.testsystems.items():
            print(testsystem)
            molecules = testsystem['molecules']

            for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator for this force field
                from simtk import openmm
                from simtk.openmm import app
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
                                                forcefield_kwargs=forcefield_kwargs,
                                                periodic_forcefield_kwargs={'nonbondedMethod':app.LJPME},
                                                nonperiodic_forcefield_kwargs={'nonbondedMethod':app.CutoffNonPeriodic},
                                                molecules=molecules)

                # Parameterize molecules
                for molecule in molecules:
                    # Create non-periodic Topology
                    nonperiodic_openmm_topology = molecule.to_topology().to_openmm()
                    system = generator.create_system(nonperiodic_openmm_topology)
                    forces = { force.__class__.__name__ : force for force in system.getForces() }
                    assert forces['NonbondedForce'].getNonbondedMethod() == openmm.NonbondedForce.CutoffNonPeriodic, "Expected CutoffNonPeriodic, got {forces['NonbondedForce'].getNonbondedMethod()}"

                    # Create periodic Topology
                    import numpy as np
                    import copy
                    box_vectors = unit.Quantity(np.diag([30, 30, 30]), unit.angstrom)
                    periodic_openmm_topology = copy.deepcopy(nonperiodic_openmm_topology)
                    periodic_openmm_topology.setPeriodicBoxVectors(box_vectors)
                    system = generator.create_system(periodic_openmm_topology)
                    forces = { force.__class__.__name__ : force for force in system.getForces() }
                    assert forces['NonbondedForce'].getNonbondedMethod() == openmm.NonbondedForce.LJPME, "Expected LJPME, got {forces['NonbondedForce'].getNonbondedMethod()}"

    def test_parameterize_molecules_from_creation(self):
        """Test that SystemGenerator can parameterize pre-specified molecules in vacuum"""
        from openmmforcefields.generators import SystemGenerator

        for name, testsystem in self.testsystems.items():
            print(testsystem)
            molecules = testsystem['molecules']

            for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator for this force field
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
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
                    assert (t2.interval() < t1.interval())

    def test_parameterize_molecules_specified_during_create_system(self):
        """Test that SystemGenerator can parameterize molecules specified during create_system"""
        from openmmforcefields.generators import SystemGenerator

        for name, testsystem in self.testsystems.items():
            molecules = testsystem['molecules']

            for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator for this force field
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield)

                # Parameterize molecules
                from openmmforcefields.utils import Timer
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    # Specify molecules during system creation
                    system = generator.create_system(openmm_topology, molecules=molecules)

    def test_add_molecules(self):
        """Test that Molecules can be added to SystemGenerator later"""
        from openmmforcefields.generators import SystemGenerator

        for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
            # Create a SystemGenerator for this force field
            generator = SystemGenerator(forcefields=self.amber_forcefields,
                                            small_molecule_forcefield=small_molecule_forcefield)

            # Add molecules for each test system separately
            for name, testsystem in self.testsystems.items():
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
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms
                    assert (t2.interval() < t1.interval())

    def test_cache(self):
        """Test that SystemGenerator correctly manages a cache"""
        from openmmforcefields.generators import SystemGenerator
        from openmmforcefields.utils import Timer

        timing = dict() # timing[(small_molecule_forcefield, smiles)] is the time (in seconds) to parameterize molecule the first time
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a single shared cache for all force fields
            cache = os.path.join(tmpdirname, 'db.json')
            # Test that we can parameterize all molecules for all test systems
            for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
                                                cache=cache)
                # Add molecules for each test system separately
                for name, testsystem in self.testsystems.items():
                    molecules = testsystem['molecules']
                    # Add molecules
                    generator.add_molecules(molecules)

                    # Parameterize molecules
                    for molecule in molecules:
                        openmm_topology = molecule.to_topology().to_openmm()
                        with Timer() as timer:
                            system = generator.create_system(openmm_topology)
                        assert system.getNumParticles() == molecule.n_atoms
                        # Record time
                        timing[(small_molecule_forcefield, molecule.to_smiles())] = timer.interval()

            # Molecules should now be cached; test timing is faster the second time
            # Test that we can parameterize all molecules for all test systems
            for small_molecule_forcefield in SystemGenerator.SMALL_MOLECULE_FORCEFIELDS:
                # Create a SystemGenerator
                generator = SystemGenerator(forcefields=self.amber_forcefields,
                                                small_molecule_forcefield=small_molecule_forcefield,
                                                molecules=molecules,
                                                cache=cache)
                # Add molecules for each test system separately
                for name, testsystem in self.testsystems.items():
                    molecules = testsystem['molecules']
                    # We don't need to add molecules that are already defined in the cache

                    # Parameterize molecules
                    for molecule in molecules:
                        openmm_topology = molecule.to_topology().to_openmm()
                        with Timer() as timer:
                            system = generator.create_system(openmm_topology)
                        assert system.getNumParticles() == molecule.n_atoms
                        # Check that this was faster than the first time
                        assert (timer.interval() < timing[(small_molecule_forcefield, molecule.to_smiles())])

    def test_complex(self):
        """Test parameterizing a protein:ligand complex in vacuum"""
        from openmmforcefields.generators import SystemGenerator
        for name, testsystem in self.testsystems.items():
            print(f'Testing parameterization of {name} in vacuum')
            molecules = testsystem['molecules']
            # Select a complex from the set
            ligand_index = 0
            complex_structure = testsystem['complex_structures'][ligand_index]
            molecule = molecules[ligand_index]
            openmm_topology = complex_structure.topology

            cache = os.path.join(get_data_filename(os.path.join('perses_jacs_systems', name)), 'cache.json')

            # Create a system in vacuum
            generator = SystemGenerator(forcefields=self.amber_forcefields,
                molecules=molecules, cache=cache)
            system = generator.create_system(openmm_topology)
            assert system.getNumParticles() == len(complex_structure.atoms)

            # Create solvated structure
            from simtk.openmm import app
            from simtk import unit
            modeller = app.Modeller(complex_structure.topology, complex_structure.positions)
            modeller.addSolvent(generator.forcefield, padding=0*unit.angstroms, ionicStrength=300*unit.millimolar)

            # Create a system with solvent and ions
            system = generator.create_system(modeller.topology)
            assert system.getNumParticles() == len(list(modeller.topology.atoms()))

            with open('test.pdb', 'w') as outfile:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
