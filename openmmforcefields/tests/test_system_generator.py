import copy
import os
import tempfile

import numpy as np
import openmm
import pytest
from openff.toolkit.topology import Molecule
from openmm.app import LJPME, PME, CutoffNonPeriodic, Modeller, PDBFile
from openmm import unit, MonteCarloBarostat, MonteCarloMembraneBarostat

from openmmforcefields.generators import SystemGenerator
from openmmforcefields.utils import Timer, get_data_filename

CI = "CI" in os.environ

################################################################################
# Tests
################################################################################


@pytest.fixture(scope="class", autouse=True)
def test_systems():
    testsystems = dict()
    for system_name, prefix in [
        # TODO: Uncomment these after we fix input files
        ("bace", "Bace"),
        # ('cdk1', 'CDK2'),
        # ('jnk1', 'Jnk1'),
        # ('mcl1', 'MCL1'),
        # ('p38', 'p38'),
        # ('ptp1b', 'PTP1B'),
        # ('thrombin', 'Thrombin'),
        # ('tyk2', 'Tyk2'),
    ]:
        # Load protein
        pdb_filename = get_data_filename(os.path.join("perses_jacs_systems", system_name, prefix + "_protein.pdb"))
        pdbfile = PDBFile(pdb_filename)

        # Load molecules
        sdf_filename = get_data_filename(
            os.path.join("perses_jacs_systems", system_name, prefix + "_ligands_shifted.sdf")
        )

        molecules = Molecule.from_file(sdf_filename, allow_undefined_stereo=True)
        print(f"Read {len(molecules)} molecules from {sdf_filename}")
        n_molecules = len(molecules)

        # Limit number of molecules for testing
        MAX_MOLECULES = 10 if not CI else 2
        if n_molecules > MAX_MOLECULES:
            print(f"Limiting to {MAX_MOLECULES} for testing...")
            n_molecules = MAX_MOLECULES
        molecules = [molecules[index] for index in range(n_molecules)]

        # Create structures
        import parmed

        # NOTE: This does not work because parmed does not correctly assign bonds for HID
        # protein_structure = parmed.load_file(pdb_filename)
        # NOTE: This is the workaround
        protein_structure = parmed.openmm.load_topology(pdbfile.topology, xyz=pdbfile.positions)

        molecules_structure = parmed.load_file(sdf_filename)
        molecules_structure = [molecules_structure[index] for index in range(n_molecules)]

        complex_structures = [(molecules_structure[index] + protein_structure) for index in range(n_molecules)]
        complex_structures = [molecules_structure[index] for index in range(n_molecules)]  # DEBUG

        # Store
        testsystem = {
            "name": system_name,
            "protein_pdbfile": pdbfile,
            "molecules": molecules,
            "complex_structures": complex_structures,
        }
        testsystems[system_name] = testsystem

        # DEBUG
        for name, testsystem in testsystems.items():
            filename = f"testsystem-{name}.pdb"
            print(filename)
            structure = testsystem["complex_structures"][0]
            # structure.save(filename, overwrite=True)
            with open(filename, "w") as outfile:
                PDBFile.writeFile(structure.topology, structure.positions, outfile)
            testsystem["molecules"][0].to_file(f"testsystem-{name}-molecule.sdf", file_format="SDF")
            testsystem["molecules"][0].to_file(f"testsystem-{name}-molecule.pdb", file_format="PDB")

    # TODO: Create other test topologies
    # TODO: Protein-only
    # TODO: Protein-ligand topology
    # TODO: Solvated protein-ligand topology
    # TODO: Host-guest topology
    # Suppress DEBUG logging from various packages

    import logging

    for name in ["parmed", "matplotlib"]:
        logging.getLogger(name).setLevel(logging.WARNING)

    return testsystems


class TestSystemGenerator:
    # AMBER force field combination to test
    amber_forcefields = [
        "amber/protein.ff14SB.xml",
        "amber/tip3p_standard.xml",
        "amber/tip3p_HFE_multivalent.xml",
    ]

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
        MAX_ATOMS = 45
        molecules = [molecule for molecule in molecules if molecule.n_atoms < MAX_ATOMS]
        # Cut down number of tests for travis
        MAX_MOLECULES = 10 if not CI else 2
        molecules = molecules[:MAX_MOLECULES]

        return molecules

    def test_create(self):
        """Test SystemGenerator creation with only OpenMM ffxml force fields"""
        # Create an empty system generator
        SystemGenerator()

    @pytest.mark.parametrize(
        "barostat_class, args",
        [
            # MonteCarloBarostat
            (
                MonteCarloBarostat,
                [0.95 * unit.atmospheres, 301.0 * unit.kelvin, 23],
            ),
            # MonteCarloMembraneBarostat
            (
                MonteCarloMembraneBarostat,
                [
                    1.0 * unit.atmospheres,
                    10.0 * unit.millinewton / unit.meter,
                    301.0 * unit.kelvin,
                    MonteCarloMembraneBarostat.XYIsotropic,
                    MonteCarloMembraneBarostat.ZFree,
                    23,
                ],
            ),
        ],
    )
    def test_barostat(self, barostat_class, args):
        """Test that different barostats are correctly applied to the system"""
        # Create a protein SystemGenerator
        generator = SystemGenerator(forcefields=self.amber_forcefields)

        # Create the barostat
        generator.barostat = barostat_class(*args)

        # Derive expected values based on barostat type
        if barostat_class is MonteCarloBarostat:
            expected = {
                "pressure": args[0],
                "temperature": args[1],
                "frequency": args[2],
            }
        else:  # MonteCarloMembraneBarostat
            expected = {
                "pressure": args[0],
                "surface_tension": args[1],
                "temperature": args[2],
                "frequency": args[-1],
            }

        # Load a PDB file
        pdb_filename = get_data_filename(os.path.join("perses_jacs_systems", "mcl1", "MCL1_protein.pdb"))
        pdbfile = PDBFile(pdb_filename)

        # Delete hydrogens from terminal protein residues
        # TODO: Fix the input files so we don't need to do this
        modeller = Modeller(pdbfile.topology, pdbfile.positions)
        residues = [residue for residue in modeller.topology.residues() if residue.name != "UNL"]
        termini_ids = [residues[0].id, residues[-1].id]

        hs = [
            atom
            for atom in modeller.topology.atoms()
            if atom.element.symbol in ["H"] and atom.residue.id in termini_ids
        ]
        modeller.delete(hs)
        modeller.addHydrogens()

        # Create a System
        system = generator.create_system(modeller.topology)

        # Check barostat is present
        forces = {force.__class__.__name__: force for force in system.getForces()}
        name = barostat_class.__name__
        assert name in forces, f"{name} not found in system forces"

        # Check barostat parameters
        force = forces[name]
        assert force.getDefaultTemperature() == expected["temperature"]
        assert force.getDefaultPressure() == expected["pressure"]
        assert force.getFrequency() == expected["frequency"]

        # Conditional check
        if hasattr(force, "getDefaultSurfaceTension"):
            assert force.getDefaultSurfaceTension() == expected["surface_tension"]

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_create_with_template_generator(self, small_molecule_forcefield):
        """Test SystemGenerator creation with small molecule residue template generators"""
        # Create a generator that defines AMBER and small molecule force fields
        generator = SystemGenerator(
            forcefields=self.amber_forcefields,
            small_molecule_forcefield=small_molecule_forcefield,
        )

        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, "db.json")
            # Create a new database file
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                cache=cache,
                small_molecule_forcefield=small_molecule_forcefield,
            )
            del generator
            # Reopen it (with cache still empty)
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                cache=cache,
                small_molecule_forcefield=small_molecule_forcefield,
            )
            del generator

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_forcefield_default_kwargs(self, small_molecule_forcefield, test_systems):
        """Test that default forcefield kwargs work correctly"""
        from openmm import unit

        forcefield_kwargs = dict()
        from openmmforcefields.generators import SystemGenerator

        for name, testsystem in test_systems.items():
            print(testsystem)
            molecules = testsystem["molecules"]

            # Create a SystemGenerator for this force field
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
                forcefield_kwargs=forcefield_kwargs,
                molecules=molecules,
            )

            # Parameterize molecules
            for molecule in molecules:
                # Create non-periodic Topology
                nonperiodic_openmm_topology = molecule.to_topology().to_openmm()
                system = generator.create_system(nonperiodic_openmm_topology)
                forces = {force.__class__.__name__: force for force in system.getForces()}
                assert forces["NonbondedForce"].getNonbondedMethod() == openmm.NonbondedForce.NoCutoff, (
                    "Expected CutoffNonPeriodic, got {forces['NonbondedForce'].getNonbondedMethod()}"
                )

                # Create periodic Topology
                box_vectors = unit.Quantity(np.diag([30, 30, 30]), unit.angstrom)
                periodic_openmm_topology = copy.deepcopy(nonperiodic_openmm_topology)
                periodic_openmm_topology.setPeriodicBoxVectors(box_vectors)
                system = generator.create_system(periodic_openmm_topology)
                forces = {force.__class__.__name__: force for force in system.getForces()}
                assert forces["NonbondedForce"].getNonbondedMethod() == openmm.NonbondedForce.PME, (
                    "Expected LJPME, got {forces['NonbondedForce'].getNonbondedMethod()}"
                )

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_forcefield_kwargs(self, small_molecule_forcefield, test_systems):
        """Test that forcefield_kwargs and nonbonded method specifications work correctly"""
        from openmm import unit

        forcefield_kwargs = {"hydrogenMass": 4 * unit.amu}
        from openmmforcefields.generators import SystemGenerator

        # Test exception is raised
        with pytest.raises(ValueError) as excinfo:
            # Not allowed to specify nonbondedMethod in forcefield_kwargs
            generator = SystemGenerator(forcefield_kwargs={"nonbondedMethod": PME})
        assert "nonbondedMethod cannot be specified in forcefield_kwargs" in str(excinfo.value)

        for name, testsystem in test_systems.items():
            print(testsystem)
            molecules = testsystem["molecules"]

            # Create a SystemGenerator for this force field
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
                forcefield_kwargs=forcefield_kwargs,
                periodic_forcefield_kwargs={"nonbondedMethod": LJPME},
                nonperiodic_forcefield_kwargs={"nonbondedMethod": CutoffNonPeriodic},
                molecules=molecules,
            )

            # Parameterize molecules
            for molecule in molecules:
                # Create non-periodic Topology
                nonperiodic_openmm_topology = molecule.to_topology().to_openmm()
                system = generator.create_system(nonperiodic_openmm_topology)
                forces = {force.__class__.__name__: force for force in system.getForces()}
                assert forces["NonbondedForce"].getNonbondedMethod() == openmm.NonbondedForce.CutoffNonPeriodic, (
                    "Expected CutoffNonPeriodic, got {forces['NonbondedForce'].getNonbondedMethod()}"
                )

                # Create periodic Topology
                box_vectors = unit.Quantity(np.diag([30, 30, 30]), unit.angstrom)
                periodic_openmm_topology = copy.deepcopy(nonperiodic_openmm_topology)
                periodic_openmm_topology.setPeriodicBoxVectors(box_vectors)
                system = generator.create_system(periodic_openmm_topology)
                forces = {force.__class__.__name__: force for force in system.getForces()}
                assert forces["NonbondedForce"].getNonbondedMethod() == openmm.NonbondedForce.LJPME, (
                    "Expected LJPME, got {forces['NonbondedForce'].getNonbondedMethod()}"
                )

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_parameterize_molecules_from_creation(self, test_systems, small_molecule_forcefield):
        """Test that SystemGenerator can parameterize pre-specified molecules in vacuum"""
        for name, testsystem in test_systems.items():
            print(testsystem)
            molecules = testsystem["molecules"]

            # Create a SystemGenerator for this force field
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
                molecules=molecules,
            )

            # Parameterize molecules
            for molecule in molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                with Timer() as t1:
                    system = generator.create_system(openmm_topology)
                assert system.getNumParticles() == molecule.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = generator.create_system(openmm_topology)
                assert system.getNumParticles() == molecule.n_atoms
                assert t2.interval() < t1.interval()

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_parameterize_molecules_specified_during_create_system(self, test_systems, small_molecule_forcefield):
        """Test that SystemGenerator can parameterize molecules specified during create_system"""
        for name, testsystem in test_systems.items():
            molecules = testsystem["molecules"]

            # Create a SystemGenerator for this force field
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
            )

            # Parameterize molecules
            for molecule in molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                # Specify molecules during system creation
                generator.create_system(openmm_topology, molecules=molecules)

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_add_molecules(self, test_systems, small_molecule_forcefield):
        """Test that Molecules can be added to SystemGenerator later"""
        # Create a SystemGenerator for this force field
        generator = SystemGenerator(
            forcefields=self.amber_forcefields,
            small_molecule_forcefield=small_molecule_forcefield,
        )

        # Add molecules for each test system separately
        for name, testsystem in test_systems.items():
            molecules = testsystem["molecules"]

            # Add molecules
            generator.add_molecules(molecules)

            # Parameterize molecules
            for molecule in molecules:
                openmm_topology = molecule.to_topology().to_openmm()
                with Timer() as t1:
                    system = generator.create_system(openmm_topology)
                assert system.getNumParticles() == molecule.n_atoms
                # Molecule should now be cached
                with Timer() as t2:
                    system = generator.create_system(openmm_topology)
                assert system.getNumParticles() == molecule.n_atoms
                assert t2.interval() < t1.interval()

    @pytest.mark.parametrize(
        "small_molecule_forcefield",
        [
            pytest.param("gaff-2.2.20", marks=pytest.mark.gaff),
            "openff-2.0.0",
            pytest.param("espaloma-0.3.2", marks=pytest.mark.espaloma),
        ],
    )
    def test_cache(self, test_systems, small_molecule_forcefield):
        """Test that SystemGenerator correctly manages a cache"""
        # timing[(small_molecule_forcefield, smiles)] is the time (in seconds) to parameterize molecule the first time
        timing = dict()
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a single shared cache for all force fields
            cache = os.path.join(tmpdirname, "db.json")
            # Test that we can parameterize all molecules for all test systems
            # Create a SystemGenerator
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
                cache=cache,
            )
            # Add molecules for each test system separately
            for name, testsystem in test_systems.items():
                molecules = testsystem["molecules"]
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
            # Create a SystemGenerator
            generator = SystemGenerator(
                forcefields=self.amber_forcefields,
                small_molecule_forcefield=small_molecule_forcefield,
                cache=cache,
            )
            # Add molecules for each test system separately
            for name, testsystem in test_systems.items():
                molecules = testsystem["molecules"]
                # We don't need to add molecules that are already defined in the cache

                # Parameterize molecules
                for molecule in molecules:
                    openmm_topology = molecule.to_topology().to_openmm()
                    with Timer() as timer:
                        system = generator.create_system(openmm_topology)
                    assert system.getNumParticles() == molecule.n_atoms

    def test_complex(self, test_systems):
        """Test parameterizing a protein:ligand complex in vacuum"""
        for name, testsystem in test_systems.items():
            from openmm import unit

            print(f"Testing parameterization of {name} in vacuum")
            molecules = testsystem["molecules"]
            # Select a complex from the set
            ligand_index = 0
            complex_structure = testsystem["complex_structures"][ligand_index]
            openmm_topology = complex_structure.topology

            cache = os.path.join(
                get_data_filename(os.path.join("perses_jacs_systems", name)),
                "cache.json",
            )

            # Create a system in vacuum
            generator = SystemGenerator(forcefields=self.amber_forcefields, molecules=molecules, cache=cache)
            system = generator.create_system(openmm_topology)
            assert system.getNumParticles() == len(complex_structure.atoms)

            # Create solvated structure
            modeller = Modeller(complex_structure.topology, complex_structure.positions)
            modeller.addSolvent(
                generator.forcefield,
                padding=0 * unit.angstroms,
                ionicStrength=300 * unit.millimolar,
            )

            # Create a system with solvent and ions
            system = generator.create_system(modeller.topology)
            assert system.getNumParticles() == len(list(modeller.topology.atoms()))

            with open("test.pdb", "w") as outfile:
                PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
