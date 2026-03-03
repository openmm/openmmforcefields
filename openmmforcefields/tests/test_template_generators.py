# ruff: noqa: E501
import glob
import random
import copy
import logging
import os
import tempfile
import unittest
import warnings
from typing import TYPE_CHECKING

import numpy as np
import openmm
import pytest
from openff.toolkit import Molecule, Topology, ForceField as OFFForceField
from openff.units import unit as OFFUnit
from openmm.app import PME, ForceField, Modeller, NoCutoff, PDBFile

from openmmforcefields.generators import (
    EspalomaTemplateGenerator,
    GAFFTemplateGenerator,
    SMIRNOFFTemplateGenerator,
)
from openmmforcefields.utils import get_data_filename

if TYPE_CHECKING:
    import parmed

_logger = logging.getLogger("openmmforcefields.tests.test_template_generators")

CI = "CI" in os.environ

CHARGE_WARNING_PATTERN = "Sum of user-provided partial charges.*does not match formal charge"

################################################################################
# Tests
################################################################################


class EnergyError(BaseException):
    pass


class ForceError(BaseException):
    pass


class TemplateGeneratorBaseCase(unittest.TestCase):
    def filter_molecules(
        self,
        molecules: list[Molecule] | list["parmed.Structure"],
        max_molecules: int = 50,
        max_atoms: int = 40,
    ) -> list[Molecule] | list["parmed.Structure"]:
        """
        Filter molecules, randomly selected from input, for the purpose of running fewer tests.

        Parameters
        ----------
        molecules
            The input list of molecules to be filtered
        max_molecules
            The maximum number of molecules to return
        max_atoms
            The maximum number of atoms in a molecule to not filter out

        Returns
        -------
        molecules
            The filtered list of molecules to be filtered
        """

        molecules = [molecule for molecule in molecules if len(molecule.atoms) <= max_atoms]

        return random.sample(molecules, min(len(molecules), max_molecules))

    def _filter_openff(self, force_field: str) -> bool:
        """
        Return True only for Parsley and force fields among SMIRNOFF force fields.

        Skips over (returns False) smirnoff99Frosst, ff14SB, water models, and 2.0.0rc1.
        """

        if "ff14sb" in force_field:
            return False
        if "tip" in force_field:
            return False
        if "opc" in force_field:
            return False
        if "spce" in force_field:
            return False

        # We cannot test openff-2.0.0-rc.1 because it triggers an openmm.OpenMMException
        # due to an equilibrium angle > \pi
        # See https://github.com/openmm/openmm/issues/3185
        if "openff" in force_field and "2.0.0-rc.1" in force_field:
            return False

        # smirnoff99Frosst is older and produces some weird geometries with some molecules
        # it's not so relevant today so skip over it in testing, though users can access it
        if "smirnoff99" in force_field:
            return False

        return True

    def setUp(self):
        from openff.units import unit

        # TODO: Harmonize with test_system_generator.py infrastructure
        # Read test molecules
        filename = get_data_filename("minidrugbank/MiniDrugBank-without-unspecified-stereochemistry.sdf")
        molecules = Molecule.from_file(filename, allow_undefined_stereo=True)

        # DEBUG: Insert acetone perturbed from planarity as first test molecule, since it fails quickly if
        # something is wrong
        molecule = Molecule.from_smiles("C=O")
        molecule.generate_conformers(n_conformers=1)

        molecule.conformers[0][0, 0] += unit.Quantity(0.1, unit.angstroms)

        molecules.insert(0, molecule)
        # DEBUG END

        # Filter molecules as appropriate
        self.molecules = self.filter_molecules(molecules, max_atoms=25 if CI else 40)

        # Get molecules for formal charge test: pick one with the lowest
        # magnitude (probably zero) and the highest magnitude (probably
        # non-zero).  Try to get different molecules if all are neutral.
        def molecule_abs_charge(molecule):
            return abs(molecule.total_charge.m_as(unit.elementary_charge))

        self.charge_test_molecules = (
            min(self.molecules, key=molecule_abs_charge),
            max(reversed(self.molecules), key=molecule_abs_charge),
        )

        # Suppress DEBUG logging from various packages
        import logging

        for name in ["parmed", "matplotlib"]:
            logging.getLogger(name).setLevel(logging.WARNING)

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
        forces = {force.__class__.__name__: force for force in system.getForces()}
        for particle_index in range(system.getNumParticles()):
            charge, sigma, epsilon = forces["NonbondedForce"].getParticleParameters(particle_index)
            system_charges.append(charge / unit.elementary_charge)
        system_charges = np.array(system_charges)

        return system_charges

    def charges_are_equal(self, system, molecule):
        """
        Return True if partial charges are equal

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

        assert system.getNumParticles() == molecule.n_atoms

        # charges_from_system returns a NumPy array that we trust to be implicitly e
        system_charges: np.ndarray = self.charges_from_system(system)

        # type(molecule.partial_charges) depends on the toolkit version
        molecule_charges: np.ndarray = molecule.partial_charges.m_as(unit.elementary_charge)

        result = np.allclose(system_charges, molecule_charges)

        if not result:
            _logger.info("Charges are not equal")
            _logger.info(f"system charges  : {system_charges}")
            _logger.info(f"molecule charges: {molecule_charges}")

        return result

    def parameterize_with_charges(self, molecule, partial_charges):
        """
        Parameterizes a molecule with given initial partial charges.

        Parameters
        ----------
        molecule : openff.toolkit.topology.Molecule
            The molecule to parameterize.  Its partial_charges will be
            overwritten.
        partial_charges
            The value to assign to the partial_charges attribute of the molecule
            before parameterization.

        Returns
        -------
        openmm.System
            The system resulting from parameterizing the molecule.
        """

        # Set up the generator.
        generator = self._make_template_generator()
        forcefield = ForceField()
        forcefield.registerTemplateGenerator(generator.generator)

        # Set the partial charges and parameterize the system.
        molecule.partial_charges = partial_charges
        generator.add_molecules(molecule)
        return forcefield.createSystem(molecule.to_topology().to_openmm(), nonbondedMethod=NoCutoff)

    def _make_template_generator(self):
        """
        Makes a template generator for generic testing.  By default, calls
        `TEMPLATE_GENERATOR` with no arguments.  Override in derived classes to
        customize this behavior.
        """

        return self.TEMPLATE_GENERATOR()

    @classmethod
    def get_permutation_indices(cls, system_1, system_2):
        """
        Computes permutations of particle indices necessary to map particles
        from one system to another, assuming that all particles and virtual
        sites are in the same order between the systems with respect to
        themselves, but not each other (i.e., particles and virtual sites may be
        interspersed among each other in different ways for the two systems).

        An exception will be raised if the numbers of particles or virtual sites
        do not match between the two systems and a permutation cannot be found.

        Parameters
        ----------
        system_1 : openmm.System
            The first system to process
        system_2 : openmm.System
            The second system to process

        Return
        ------
        (list[int], list[int])
            A permutation (for each particle in the first system, the index of
            the particle in the second system) and its inverse.
        """

        atoms_1, sites_1, atoms_2, sites_2 = [], [], [], []
        for index in range(system_1.getNumParticles()):
            if system_1.isVirtualSite(index):
                sites_1.append(index)
            else:
                atoms_1.append(index)
        for index in range(system_2.getNumParticles()):
            if system_2.isVirtualSite(index):
                sites_2.append(index)
            else:
                atoms_2.append(index)

        assert len(atoms_1) == len(atoms_2)
        assert len(sites_1) == len(sites_2)

        particles_1 = atoms_1 + sites_1
        particles_2 = atoms_2 + sites_2
        permutation_12 = [-1] * len(particles_1)
        permutation_21 = [-1] * len(particles_1)
        for index_1, index_2 in zip(particles_1, particles_2):
            permutation_12[index_1] = index_2
            permutation_21[index_2] = index_1
        return permutation_12, permutation_21

    @classmethod
    def compare_energies(
        cls,
        molecule_name: str,
        positions: openmm.unit.Quantity,
        template_generated_system: openmm.System,
        reference_system: openmm.System,
        extra_info: str = "",
    ):
        """
        Compare energies between OpenMM System generated by reference method and OpenMM System generated
        by ForceField template.

        The OpenMM System object created internally by the reference method is used to
        avoid any issues with stochasticity of partial charges due to conformer generation.

        Parameters
        ----------
        molecule_name : str
            A descriptive name for the molecule (can be its SMILES string) to put in error messages
        positions : openmm.unit.Quantity or None
            Positions for the particles in the template generated system
        template_generated_system : openmm.System
            System generated by OpenMM ForceField template
        reference_system : openmm.System
            System generated by reference parmaeterization engine
        extra_info : str
            Extra information to include in the error message
        """

        ENERGY_UNIT = openmm.unit.kilocalories_per_mole
        FORCE_UNIT = openmm.unit.kilocalories_per_mole / openmm.unit.angstroms

        # Get permutations to handle differing orders of virtual sites
        forces_permutation, positions_permutation = cls.get_permutation_indices(
            template_generated_system, reference_system
        )

        # Sanity check: make sure that the correct number of constraints are present
        assert template_generated_system.getNumConstraints() == reference_system.getNumConstraints()

        # Compute energies and forces
        template_energy, template_forces = cls.compute_energy(template_generated_system, positions)
        reference_positions = [positions[index] for index in positions_permutation]
        reference_energy, reference_forces = cls.compute_energy(reference_system, reference_positions)
        forces_total = reference_forces["total"]
        reference_forces["total"] = (
            np.array([forces_total[index].value_in_unit(FORCE_UNIT) for index in forces_permutation]) * FORCE_UNIT
        )
        for component, forces_component in reference_forces["components"].items():
            reference_forces["components"][component] = (
                np.array([forces_component[index].value_in_unit(FORCE_UNIT) for index in forces_permutation])
                * FORCE_UNIT
            )

        def write_xml(filename, system):
            with open(filename, "w") as outfile:
                print(f"Writing {filename}...")
                outfile.write(openmm.XmlSerializer.serialize(system))

        # Make sure both systems contain the same energy components
        reference_components = set(reference_energy["components"])
        template_components = set(template_energy["components"])
        if len(reference_components.difference(template_components)) > 0:
            raise Exception(
                f"Reference system contains components {reference_components.difference(template_components)} "
                "that do not appear in template-generated system."
            )
        if len(template_components.difference(reference_components)) > 0:
            raise Exception(
                "Template-generated system contains components "
                f"{template_components.difference(reference_components)} "
                "that do not appear in reference system."
            )
        components = reference_components

        # Compare energies
        ENERGY_DEVIATION_TOLERANCE = 1.0e-2 * ENERGY_UNIT
        delta = template_energy["total"] - reference_energy["total"]
        if abs(delta) > ENERGY_DEVIATION_TOLERANCE:
            # Show breakdown by components
            print("Energy components:")
            print(f"{'component':24} {'Template (kcal/mol)':>20} {'Reference (kcal/mol)':>20}")
            for key in components:
                reference_component_energy = reference_energy["components"][key].value_in_unit(ENERGY_UNIT)
                template_component_energy = template_energy["components"][key].value_in_unit(ENERGY_UNIT)
                print(f"{key:24} {template_component_energy:20.3f} {reference_component_energy:20.3f}")

            reference_total_energy = reference_energy["total"].value_in_unit(ENERGY_UNIT)
            template_total_energy = template_energy["total"].value_in_unit(ENERGY_UNIT)
            print(f"{'TOTAL':24} {template_total_energy:20.3f} {reference_total_energy:20.3f}")
            write_xml("reference_system.xml", reference_system)
            write_xml("test_system.xml", template_generated_system)
            raise EnergyError(
                f"Energy deviation for {molecule_name} ({delta.in_units_of(ENERGY_UNIT)}) exceeds "
                f"threshold ({ENERGY_DEVIATION_TOLERANCE}). {extra_info=}"
            )

        # Compare forces
        def norm_sq(x):
            return (x * x).sum() / x.shape[0]

        def relative_deviation(x, y):
            x = x.value_in_unit(FORCE_UNIT)
            y = y.value_in_unit(FORCE_UNIT)
            scale = norm_sq(x) + norm_sq(y)
            return np.sqrt(norm_sq(x - y) / scale) if scale > 0 else 0

        RELATIVE_FORCE_DEVIATION_TOLERANCE = 1.0e-5
        relative_force_deviation = relative_deviation(template_forces["total"], reference_forces["total"])
        if relative_force_deviation > RELATIVE_FORCE_DEVIATION_TOLERANCE:
            # Show breakdown by components
            print("Force components:")
            print(f"{'component':24} {'relative deviation':>24}")
            for key in components:
                deviation = relative_deviation(template_forces["components"][key], reference_forces["components"][key])
                print(f"{key:24} {deviation:24.10f}")

            print(f"{'TOTAL':24} {relative_force_deviation:24.10f}")
            write_xml("reference_system.xml", reference_system)
            write_xml("test_system.xml", template_generated_system)
            raise ForceError(
                f"Relative force deviation for {molecule_name} ({relative_force_deviation}) exceeds threshold "
                f"({RELATIVE_FORCE_DEVIATION_TOLERANCE}). {extra_info=}"
            )

    @classmethod
    def compare_energies_single(
        cls,
        molecule: Molecule,
        template_generated_system: openmm.System,
        reference_system: openmm.System,
        extra_info: str = "",
    ):
        """
        Compare energies between OpenMM System generated by reference method and OpenMM System generated
        by ForceField template.  Wrapper around `compare_energies` for the special case of a single molecule
        with no virtual sites.  The test positions are taken from conformer 0 of the molecule provided.
        """

        return cls.compare_energies(
            molecule.to_smiles(),
            molecule.conformers[0].to_openmm(),
            template_generated_system,
            reference_system,
            extra_info,
        )

    @staticmethod
    def compute_energy(system, positions):
        """
        Compute potential energy and Force components for an OpenMM system.

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
        platform = openmm.Platform.getPlatformByName("Reference")
        integrator = openmm.VerletIntegrator(0.001)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        context.applyConstraints(integrator.getConstraintTolerance())

        openmm_energy = {
            "total": context.getState(getEnergy=True).getPotentialEnergy(),
            "components": {
                system.getForce(index).__class__.__name__: context.getState(
                    getEnergy=True, groups=(1 << index)
                ).getPotentialEnergy()
                for index in range(system.getNumForces())
            },
        }

        openmm_forces = {
            "total": context.getState(getForces=True).getForces(asNumpy=True),
            "components": {
                system.getForce(index).__class__.__name__: context.getState(
                    getForces=True, groups=(1 << index)
                ).getForces(asNumpy=True)
                for index in range(system.getNumForces())
            },
        }

        del context, integrator
        return openmm_energy, openmm_forces

    def propagate_dynamics(self, positions, system):
        """
        Run a few steps of dynamics to generate a perturbed configuration.

        Parameters
        ----------
        positions : openmm.unit.Quantity
            Initial positions to start dynamics from
        system : openmm.System
            System object for dynamics

        Returns
        -------
        new_positions: openmm.unit.Quantity
            Updated positions after dynamics
        """

        import openmm.unit

        temperature = 300 * openmm.unit.kelvin
        collision_rate = 1.0 / openmm.unit.picoseconds
        timestep = 1.0 * openmm.unit.femtoseconds
        nsteps = 100
        integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
        platform = openmm.Platform.getPlatformByName("Reference")
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)

        # RDKit-generated conformer is occassionally a bad starting point for dynamics,
        # so minimize with MM first, see https://github.com/openmm/openmmforcefields/pull/370#issuecomment-2749237209
        openmm.LocalEnergyMinimizer.minimize(context)
        integrator.step(nsteps)

        return context.getState(getPositions=True).getPositions()

    def propagate_dynamics_single(self, molecule, system):
        """
        Calls `propagate_dynamics` using a molecule containing positions.

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

        from openff.units.openmm import from_openmm

        # Run some dynamics
        new_positions = self.propagate_dynamics(molecule.conformers[0].to_openmm(), system)

        # Copy the molecule, storing new conformer
        new_molecule = copy.deepcopy(molecule)
        new_molecule.conformers[0] = from_openmm(new_positions)

        return new_molecule


@pytest.mark.gaff
class TestGAFFTemplateGenerator(TemplateGeneratorBaseCase):
    TEMPLATE_GENERATOR = GAFFTemplateGenerator

    amber_forcefields = [
        "amber/protein.ff14SB.xml",
        "amber/tip3p_standard.xml",
        "amber/tip3p_HFE_multivalent.xml",
    ]

    def test_version(self):
        """Test version"""
        for forcefield in GAFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator = GAFFTemplateGenerator(forcefield=forcefield)
            import re

            result = re.match(r"^gaff-(?P<major_version>\d+)\.(?P<minor_version>[\d.]+)$", forcefield)
            assert generator.forcefield == forcefield
            assert generator.gaff_version == result["major_version"] + "." + result["minor_version"]
            assert generator.gaff_major_version == result["major_version"]
            assert generator.gaff_minor_version == result["minor_version"]
            assert generator.gaff_dat_filename.endswith(forcefield + ".dat")
            assert os.path.exists(generator.gaff_dat_filename)
            assert generator.gaff_xml_filename.endswith(forcefield + ".xml")
            assert os.path.exists(generator.gaff_xml_filename)

    def test_create(self):
        """Test template generator creation"""
        # Create an empty generator
        generator = self.TEMPLATE_GENERATOR()
        # Create a generator that knows about a few molecules
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        # Create a generator that also has a database cache
        with tempfile.TemporaryDirectory() as tmpdirname:
            cache = os.path.join(tmpdirname, "db.json")
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
            assert str(e).startswith("No template found")

        # Now add the molecule to the generator and ensure parameterization passes
        generator.add_molecules(molecule)
        openmm_topology = molecule.to_topology().to_openmm()
        try:
            system = forcefield.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
        except Exception as e:
            print(forcefield._atomTypes.keys())

            PDBFile.writeFile(
                openmm_topology,
                molecule.conformers[0].to_openmm(),
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

    def test_charge_none(self):
        """Test that charges are nonzero after charging if the molecule has None for user charges"""
        for molecule in self.charge_test_molecules:
            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(molecule, None)

            # Molecule should have partial_charges set to None but system should
            # have non-zero charges assigned
            assert molecule.partial_charges is None
            assert not np.allclose(self.charges_from_system(system), 0)

    def test_charge_zero(self):
        """Test that charges are nonzero after charging if the molecule has zero for user charges"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(
                    molecule, unit.Quantity(np.zeros(molecule.n_atoms), unit.elementary_charge)
                )

            # Molecule should have partial_charges set to zero but system should
            # have non-zero charges assigned
            assert molecule.partial_charges is not None
            assert np.all(molecule.partial_charges.m_as(unit.elementary_charge) == 0)
            assert not np.allclose(self.charges_from_system(system), 0)

    def test_charge_valid(self):
        """Test that user-specified partial charges are used if requested"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            user_charges = (
                np.linspace(-1, 1, molecule.n_atoms)
                + molecule.total_charge.m_as(unit.elementary_charge) / molecule.n_atoms
            )

            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(molecule, unit.Quantity(user_charges, unit.elementary_charge))

            # User, molecule, and system charges should be equal
            assert molecule.partial_charges is not None
            assert self.charges_are_equal(system, molecule)
            assert np.allclose(self.charges_from_system(system), user_charges)

    def test_charge_warning(self):
        """Test that a warning is raised with user-specified partial charges that do not sum to the total charge"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            user_charges = (
                np.linspace(-1, 1, molecule.n_atoms)
                + molecule.total_charge.m_as(unit.elementary_charge) / molecule.n_atoms
                + 1
            )

            # Charge mismatch warning should be raised
            with pytest.warns(match=CHARGE_WARNING_PATTERN):
                system = self.parameterize_with_charges(molecule, unit.Quantity(user_charges, unit.elementary_charge))

            # User, molecule, and system charges should be equal
            assert molecule.partial_charges is not None
            assert self.charges_are_equal(system, molecule)
            assert np.allclose(self.charges_from_system(system), user_charges)

    def test_debug_ffxml(self):
        """Test that debug ffxml file is created when requested"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            debug_ffxml_filename = os.path.join(tmpdirname, "molecule.ffxml")
            cache = os.path.join(tmpdirname, "db.json")
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
            if hasattr(generator, "gaff_xml_filename"):
                forcefield_from_ffxml.loadFile(generator.gaff_xml_filename)
            forcefield_from_ffxml.loadFile(debug_ffxml_filename)
            system2 = forcefield_from_ffxml.createSystem(openmm_topology, nonbondedMethod=NoCutoff)
            # TODO: Test that systems are equivalent
            assert system.getNumParticles() == system2.getNumParticles()

    def test_cache(self):
        """Test template generator cache capability"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Create a generator that also has a database cache
            cache = os.path.join(tmpdirname, "db.json")
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
                assert n_entries == n_expected, (
                    f"Expected {n_expected} entries but database has {n_entries}\n db contents: {db_entries}"
                )

            check_cache(generator, len(self.molecules))

            # Clean up, forcing closure of database
            del forcefield, generator

            # Create a generator that also uses the database cache but has no molecules
            print("Creating new generator with just cache...")
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
        from openmm import unit

        molecule = self.molecules[0]
        openmm_topology = molecule.to_topology().to_openmm()
        openmm_positions = molecule.conformers[0].to_openmm()

        # Try adding solvent without residue template generator; this will fail
        forcefield = ForceField("tip3p.xml")
        # Add solvent to a system containing a small molecule
        modeller = Modeller(openmm_topology, openmm_positions)
        try:
            modeller.addSolvent(forcefield, model="tip3p", padding=6.0 * unit.angstroms)
        except ValueError:
            pass

        # Create a generator that knows about a few molecules
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        # Add to the forcefield object
        forcefield.registerTemplateGenerator(generator.generator)
        # Add solvent to a system containing a small molecule
        # This should succeed
        modeller.addSolvent(forcefield, model="tip3p", padding=6.0 * unit.angstroms)

    def test_jacs_ligands(self):
        """Use template generator to parameterize the Schrodinger JACS set of ligands"""
        jacs_systems = {
            # 'bace'     : { 'prefix' : 'Bace' },
            # 'cdk2'     : { 'prefix' : 'CDK2' },
            "jnk1": {"prefix": "Jnk1"},
            "mcl1": {"prefix": "MCL1"},
            # 'p38'      : { 'prefix' : 'p38' },
            "ptp1b": {"prefix": "PTP1B"},
            "thrombin": {"prefix": "Thrombin"},
            # 'tyk2'     : { 'prefix' : 'Tyk2' },
        }
        for system_name in jacs_systems:
            prefix = jacs_systems[system_name]["prefix"]
            # Load molecules
            ligand_sdf_filename = get_data_filename(
                os.path.join("perses_jacs_systems", system_name, prefix + "_ligands.sdf")
            )
            print(f"Reading molecules from {ligand_sdf_filename} ...")
            molecules = Molecule.from_file(ligand_sdf_filename, allow_undefined_stereo=True)
            # Ensure this is a list
            try:
                len(molecules)
            except TypeError:
                molecules = [molecules]

            print(f"Read {len(molecules)} molecules from {ligand_sdf_filename}")

            molecules = self.filter_molecules(molecules, max_molecules=3 if CI else len(molecules))

            print(f"{len(molecules)} molecules remain after filtering")

            # Create template generator with local cache
            cache = os.path.join(
                get_data_filename(os.path.join("perses_jacs_systems", system_name)),
                "cache.json",
            )
            generator = self.TEMPLATE_GENERATOR(molecules=molecules, cache=cache)

            # Create a ForceField
            forcefield = ForceField()
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)

            # Parameterize all molecules
            print(f"Caching all molecules for {system_name} at {cache} ...")
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
            print(f"{n_failure}/{n_success + n_failure} ligands failed to parameterize for {system_name}")

    def test_jacs_complexes(self):
        """Use template generator to parameterize the Schrodinger JACS set of complexes"""
        # TODO: Uncomment working systems when we have cleaned up the input files
        jacs_systems = {
            # 'bace'     : { 'prefix' : 'Bace' },
            # 'cdk2'     : { 'prefix' : 'CDK2' },
            # 'jnk1'     : { 'prefix' : 'Jnk1' },
            "mcl1": {"prefix": "MCL1"},
            # 'p38'      : { 'prefix' : 'p38' },
            # 'ptp1b'    : { 'prefix' : 'PTP1B' },
            # 'thrombin' : { 'prefix' : 'Thrombin' },
            # 'tyk2'     : { 'prefix' : 'Tyk2' },
        }
        for system_name in jacs_systems:
            prefix = jacs_systems[system_name]["prefix"]
            # Read molecules
            ligand_sdf_filename = get_data_filename(
                os.path.join("perses_jacs_systems", system_name, prefix + "_ligands.sdf")
            )
            print(f"Reading molecules from {ligand_sdf_filename} ...")
            molecules = Molecule.from_file(ligand_sdf_filename, allow_undefined_stereo=True)
            try:
                len(molecules)
            except TypeError:
                molecules = [molecules]
            print(f"Read {len(molecules)} molecules from {ligand_sdf_filename}")

            # Read ParmEd Structures
            import parmed
            from openmm import unit

            protein_pdb_filename = get_data_filename(
                os.path.join("perses_jacs_systems", system_name, prefix + "_protein.pdb")
            )
            print(f"Reading protein from {protein_pdb_filename} ...")

            pdbfile = PDBFile(protein_pdb_filename)
            protein_structure = parmed.openmm.load_topology(
                pdbfile.topology, xyz=pdbfile.positions.value_in_unit(unit.angstroms)
            )
            ligand_structures = parmed.load_file(ligand_sdf_filename)
            try:
                len(ligand_structures)
            except TypeError:
                ligand_structures = [ligand_structures]
            assert len(ligand_structures) == len(molecules)

            ligand_structures = self.filter_molecules(ligand_structures, max_molecules=3 if CI else 6, max_atoms=50)

            print(f"{len(ligand_structures)} molecules remain after filtering")

            # Create complexes
            complex_structures = [(protein_structure + ligand_structure) for ligand_structure in ligand_structures]

            # Create template generator with local cache
            cache = os.path.join(
                get_data_filename(os.path.join("perses_jacs_systems", system_name)),
                "cache.json",
            )
            generator = self.TEMPLATE_GENERATOR(molecules=molecules, cache=cache)

            # Create a ForceField
            forcefield = ForceField(*self.amber_forcefields)
            # Register the template generator
            forcefield.registerTemplateGenerator(generator.generator)

            # Parameterize all complexes
            print(f"Caching all molecules for {system_name} at {cache} ...")
            for ligand_index, complex_structure in enumerate(complex_structures):
                molecule = molecules[ligand_index]

                # Delete hydrogens from terminal protein residues
                # TODO: Fix the input files so we don't need to do this
                from openmm import app

                modeller = app.Modeller(complex_structure.topology, complex_structure.positions)
                residues = [residue for residue in modeller.topology.residues() if residue.name != "UNL"]
                termini_ids = [residues[0].id, residues[-1].id]

                hs = [
                    atom
                    for atom in modeller.topology.atoms()
                    if atom.element.symbol in ["H"] and atom.residue.id in termini_ids
                ]
                modeller.delete(hs)
                modeller.addHydrogens(forcefield)

                # Parameterize protein:ligand complex in vacuum
                print(f" Parameterizing {system_name} : {molecule.to_smiles()} in vacuum...")
                forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

                # Parameterize protein:ligand complex in solvent
                print(f" Parameterizing {system_name} : {molecule.to_smiles()} in explicit solvent...")
                modeller.addSolvent(
                    forcefield,
                    padding=0 * unit.angstroms,
                    ionicStrength=300 * unit.millimolar,
                )
                forcefield.createSystem(modeller.topology, nonbondedMethod=PME)

    def test_parameterize(self):
        """Test parameterizing molecules with template generator for all supported force fields"""
        # Test all supported small molecule force fields
        for small_molecule_forcefield in self.TEMPLATE_GENERATOR.INSTALLED_FORCEFIELDS:
            if not self._filter_openff(small_molecule_forcefield):
                _logger.debug(f"skipping {small_molecule_forcefield}")
                continue
            if "spce" in small_molecule_forcefield:
                continue

            _logger.info(f"Testing {small_molecule_forcefield}")
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
                assert t2.interval() < t1.interval()

    def test_multiple_registration(self):
        """Test registering the template generator with multiple force fields"""
        generator = self.TEMPLATE_GENERATOR(molecules=self.molecules)
        NUM_FORCEFIELDS = 2  # number of force fields to test
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


class TestSMIRNOFFTemplateGenerator(TemplateGeneratorBaseCase):
    TEMPLATE_GENERATOR = SMIRNOFFTemplateGenerator

    def _make_template_generator(self):
        """
        Makes a `SMIRNOFFTemplateGenerator` for generic testing.
        """

        return SMIRNOFFTemplateGenerator(forcefield="openff-2.3.0")

    def test_INSTALLED_FORCEFIELDS(self):
        """Test that names in INSTALLED_FORCEFIELDS resolve correctly"""

        assert sorted(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS) == sorted(
            SMIRNOFFTemplateGenerator._INSTALLED_FORCEFIELDS
        )

        for name, path in SMIRNOFFTemplateGenerator._INSTALLED_FORCEFIELDS.items():
            generator_name = SMIRNOFFTemplateGenerator(forcefield=name)
            generator_path = SMIRNOFFTemplateGenerator(forcefield=path)

            assert generator_name.smirnoff_filenames == generator_path.smirnoff_filenames

            assert generator_name.forcefield == name
            assert len(generator_name.smirnoff_filenames) == 1
            assert generator_name.smirnoff_filenames[0].endswith(".offxml")
            assert os.path.exists(generator_name.smirnoff_filenames[0])

    def test_forcefield_no_default(self):
        """Test that not specifying a force field gives an error"""

        with pytest.raises(ValueError):
            SMIRNOFFTemplateGenerator()

    def test_forcefield_path(self):
        """Test that specifying a path to a force field loads that force field"""

        for forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator_ref = SMIRNOFFTemplateGenerator(forcefield=forcefield)
            generator_test = SMIRNOFFTemplateGenerator(forcefield=generator_ref.smirnoff_filenames[0])
            assert generator_test.forcefield
            assert generator_test.smirnoff_filenames == generator_ref.smirnoff_filenames
            assert generator_test._smirnoff_forcefield.to_string() == generator_ref._smirnoff_forcefield.to_string()

    def test_forcefield_file(self):
        """Test that a force field can be loaded directly from a file object"""

        for forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator_ref = SMIRNOFFTemplateGenerator(forcefield=forcefield)
            with open(generator_ref.smirnoff_filenames[0]) as file:
                generator_test = SMIRNOFFTemplateGenerator(forcefield=file)
            assert generator_test.forcefield
            assert generator_test.smirnoff_filenames == [None]
            assert generator_ref._smirnoff_forcefield.to_string() == generator_test._smirnoff_forcefield.to_string()

    def test_forcefield_str(self):
        """Test that a force field can be parsed directly from a string"""

        for forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator_ref = SMIRNOFFTemplateGenerator(forcefield=forcefield)
            with open(generator_ref.smirnoff_filenames[0]) as file:
                generator_test = SMIRNOFFTemplateGenerator(forcefield=file.read())
            assert generator_test.forcefield
            assert generator_test.smirnoff_filenames == [None]
            assert generator_ref._smirnoff_forcefield.to_string() == generator_test._smirnoff_forcefield.to_string()

    def test_forcefield_bytes(self):
        """Test that a force field can be parsed directly from bytes"""

        for forcefield in SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS:
            generator_ref = SMIRNOFFTemplateGenerator(forcefield=forcefield)
            with open(generator_ref.smirnoff_filenames[0], "rb") as file:
                generator_test = SMIRNOFFTemplateGenerator(forcefield=file.read())
            assert generator_test.forcefield
            assert generator_test.smirnoff_filenames == [None]
            assert generator_ref._smirnoff_forcefield.to_string() == generator_test._smirnoff_forcefield.to_string()

    def test_forcefield_multiple(self):
        """Test loading multiple force field components together"""

        generator = SMIRNOFFTemplateGenerator(
            forcefield=["openff-2.0.0.offxml", "tip3p.offxml", OFFForceField().to_string()]
        )
        assert generator.forcefield
        assert len(generator.smirnoff_filenames) == 3
        assert generator.smirnoff_filenames[0].endswith("openff-2.0.0.offxml")
        assert generator.smirnoff_filenames[1].endswith("tip3p.offxml")
        assert generator.smirnoff_filenames[2] is None
        assert os.path.exists(generator.smirnoff_filenames[0])
        assert os.path.exists(generator.smirnoff_filenames[1])

    def test_forcefield_unconstrained(self):
        """Test that forcefield names resolve to unconstrained versions correctly"""

        # Should redirect to the unconstrained version
        generator = SMIRNOFFTemplateGenerator(forcefield="openff-2.0.0")
        assert generator.forcefield == "openff-2.0.0"
        assert len(generator.smirnoff_filenames) == 1
        assert generator.smirnoff_filenames[0].endswith("openff_unconstrained-2.0.0.offxml")
        assert os.path.exists(generator.smirnoff_filenames[0])

        # Should not point to any valid force field
        with pytest.raises(ValueError, match="Can't load or parse specified SMIRNOFF"):
            SMIRNOFFTemplateGenerator(forcefield="openff_unconstrained-2.0.0")

        # Specifying a name with .offxml should never get modified
        generator = SMIRNOFFTemplateGenerator(forcefield="openff-2.0.0.offxml")
        assert generator.forcefield
        assert len(generator.smirnoff_filenames) == 1
        assert generator.smirnoff_filenames[0].endswith("openff-2.0.0.offxml")
        assert os.path.exists(generator.smirnoff_filenames[0])

        generator = SMIRNOFFTemplateGenerator(forcefield="openff_unconstrained-2.0.0.offxml")
        assert generator.forcefield
        assert len(generator.smirnoff_filenames) == 1
        assert generator.smirnoff_filenames[0].endswith("openff_unconstrained-2.0.0.offxml")
        assert os.path.exists(generator.smirnoff_filenames[0])

    def test_energies(self):
        """Test potential energies match between openff-toolkit and OpenMM ForceField"""

        from openff.toolkit import get_available_force_fields

        # Test all SMIRNOFF force fields
        for small_molecule_forcefield in get_available_force_fields():
            if not self._filter_openff(small_molecule_forcefield):
                _logger.debug(f"skipping {small_molecule_forcefield}")
                continue
            if small_molecule_forcefield == "openff_unconstrained-2.3.0.offxml":
                # OpenFF FF with NAGL: test all molecules in the set
                molecules = self.molecules
            elif small_molecule_forcefield in SMIRNOFFTemplateGenerator._INSTALLED_FORCEFIELDS.values():
                # Other preset force field: test a subset of molecules
                molecules = self.filter_molecules(self.molecules, max_molecules=5)
            else:
                # Something else (pre-release, etc.): just try one molecule to ensure it doesn't fail
                molecules = [Molecule.from_smiles("C=O")]
                molecules[0].generate_conformers(n_conformers=1)

            _logger.info(f"Testing {small_molecule_forcefield}")
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Parameterize some molecules
            for molecule in molecules:
                # Create OpenMM System using OpenMM app
                openmm_system = openmm_forcefield.createSystem(
                    molecule.to_topology().to_openmm(),
                    removeCMMotion=False,
                    nonbondedMethod=NoCutoff,
                )

                # Retrieve System generated by the SMIRNOFF typing engine
                smirnoff_system = generator.get_openmm_system(molecule)

                # Compare energies and forces
                self.compare_energies_single(
                    molecule, openmm_system, smirnoff_system, f"uses {small_molecule_forcefield}"
                )

                # Run some dynamics
                molecule = self.propagate_dynamics_single(molecule, smirnoff_system)

                # Compare energies again
                self.compare_energies_single(
                    molecule, openmm_system, smirnoff_system, f"uses {small_molecule_forcefield}"
                )

    def test_partial_charges_are_none(self):
        """Test parameterizing a small molecule with `partial_charges=None` instead
        of zeros (happens frequently in OFFTK>=0.7.0)"""
        from openff.units import unit
        from openff.toolkit import get_available_force_fields

        molecule = Molecule.from_smiles("C=O")
        molecule.generate_conformers(n_conformers=1)
        # molecule._partial_charges = None
        assert (molecule.partial_charges is None) or np.all(molecule.partial_charges / unit.elementary_charge == 0)
        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in get_available_force_fields():
            if not self._filter_openff(small_molecule_forcefield):
                _logger.debug(f"skipping {small_molecule_forcefield}")
                continue

            _logger.info(f"Testing {small_molecule_forcefield}")
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = SMIRNOFFTemplateGenerator(molecules=[molecule], forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Create OpenMM System using OpenMM app
            openmm_forcefield.createSystem(
                molecule.to_topology().to_openmm(),
                removeCMMotion=False,
                nonbondedMethod=NoCutoff,
            )
            generator.get_openmm_system(molecule)

    def get_terms(self, system):
        """
        Helper for `test_constraints()` and `test_constraints_water()`.  Returns
        the set of harmonic bonds, harmonic angles, and constraints from a
        system.
        """

        bonds = set()
        angles = set()
        constraints = set()

        for force in system.getForces():
            if isinstance(force, openmm.HarmonicBondForce):
                for i_bond in range(force.getNumBonds()):
                    i, j = force.getBondParameters(i_bond)[:2]
                    bonds.add(min((i, j), (j, i)))

            if isinstance(force, openmm.HarmonicAngleForce):
                for i_angle in range(force.getNumAngles()):
                    i, j, k = force.getAngleParameters(i_angle)[:3]
                    angles.add(min((i, j, k), (k, j, i)))

        for i_constraint in range(system.getNumConstraints()):
            i, j = system.getConstraintParameters(i_constraint)[:2]
            constraints.add(min((i, j), (j, i)))

        return bonds, angles, constraints

    def test_constraints(self):
        """
        Tests that the expected constraints are added using the unconstrained
        and constrained variants of a SMIRNOFF force field, and OpenMM's
        `constraints` and `flexibleConstraints` options.
        """

        from openmm.app import HBonds, AllBonds, HAngles

        # Make acetaldehyde (mapped SMILES guarantees atom order)
        molecule = Molecule.from_mapped_smiles("[C:1]([C:2](=[O:3])[H:7])([H:4])([H:5])[H:6]")
        topology = molecule.to_topology().to_openmm()

        # Expected bonds, angles, and constraints:
        non_h_bonds = {(0, 1), (1, 2)}
        h_bonds = {(0, 3), (0, 4), (0, 5), (1, 6)}
        non_h_angles = {(0, 1, 2), (0, 1, 6), (1, 0, 3), (1, 0, 4), (1, 0, 5), (2, 1, 6)}
        h_angles = {(3, 0, 4), (3, 0, 5), (4, 0, 5)}
        h_angles_constraints = {(3, 4), (3, 5), (4, 5)}

        # Make force fields
        unconstrained = openmm.app.ForceField()
        unconstrained.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff_unconstrained-2.3.0.offxml").generator
        )
        constrained = openmm.app.ForceField()
        constrained.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff-2.3.0.offxml").generator
        )

        # Unconstrained force field
        assert self.get_terms(unconstrained.createSystem(topology, constraints=None)) == (
            non_h_bonds | h_bonds,
            non_h_angles | h_angles,
            set(),
        )
        assert self.get_terms(unconstrained.createSystem(topology, constraints=HBonds)) == (
            non_h_bonds,
            non_h_angles | h_angles,
            h_bonds,
        )
        assert self.get_terms(unconstrained.createSystem(topology, constraints=AllBonds)) == (
            set(),
            non_h_angles | h_angles,
            non_h_bonds | h_bonds,
        )
        assert self.get_terms(unconstrained.createSystem(topology, constraints=HAngles)) == (
            set(),
            non_h_angles,
            non_h_bonds | h_bonds | h_angles_constraints,
        )

        # Unconstrained force field with flexibleConstraints: should still get harmonic terms
        assert self.get_terms(unconstrained.createSystem(topology, constraints=HBonds, flexibleConstraints=True)) == (
            non_h_bonds | h_bonds,
            non_h_angles | h_angles,
            h_bonds,
        )
        assert self.get_terms(
            unconstrained.createSystem(topology, constraints=AllBonds, flexibleConstraints=True)
        ) == (non_h_bonds | h_bonds, non_h_angles | h_angles, non_h_bonds | h_bonds)
        assert self.get_terms(unconstrained.createSystem(topology, constraints=HAngles, flexibleConstraints=True)) == (
            non_h_bonds | h_bonds,
            non_h_angles | h_angles,
            non_h_bonds | h_bonds | h_angles_constraints,
        )

        # Constrained force field: constraints up to HBonds should be forced on
        assert self.get_terms(constrained.createSystem(topology, constraints=None)) == (
            non_h_bonds,
            non_h_angles | h_angles,
            h_bonds,
        )
        assert self.get_terms(constrained.createSystem(topology, constraints=HBonds)) == (
            non_h_bonds,
            non_h_angles | h_angles,
            h_bonds,
        )
        assert self.get_terms(constrained.createSystem(topology, constraints=AllBonds)) == (
            set(),
            non_h_angles | h_angles,
            non_h_bonds | h_bonds,
        )
        assert self.get_terms(constrained.createSystem(topology, constraints=HAngles)) == (
            set(),
            non_h_angles,
            non_h_bonds | h_bonds | h_angles_constraints,
        )

        # Constrained force field with flexibleConstraints: should not have harmonic terms for enforced constraints
        assert self.get_terms(constrained.createSystem(topology, constraints=AllBonds, flexibleConstraints=True)) == (
            non_h_bonds,
            non_h_angles | h_angles,
            non_h_bonds | h_bonds,
        )
        assert self.get_terms(constrained.createSystem(topology, constraints=HAngles, flexibleConstraints=True)) == (
            non_h_bonds,
            non_h_angles | h_angles,
            non_h_bonds | h_bonds | h_angles_constraints,
        )

    def test_unconstrained_default(self):
        """
        Tests that force field names by default select unconstrained variants.
        """

        molecule = Molecule.from_mapped_smiles("[C:1]([C:2](=[O:3])[H:7])([H:4])([H:5])[H:6]")
        topology = molecule.to_topology().to_openmm()

        unconstrained = openmm.app.ForceField()
        unconstrained.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff_unconstrained-2.3.0.offxml").generator
        )
        constrained = openmm.app.ForceField()
        constrained.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff-2.3.0.offxml").generator
        )
        default = openmm.app.ForceField()
        default.registerTemplateGenerator(
            SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff-2.3.0").generator
        )

        assert unconstrained.createSystem(topology).getNumConstraints() == 0
        assert constrained.createSystem(topology).getNumConstraints() > 0
        assert default.createSystem(topology).getNumConstraints() == 0

    def test_constraints_water(self):
        """
        Tests that the expected constraints are added for SMIRNOFF water using
        OpenMM's `rigidWater` option.
        """

        # Make water (mapped SMILES guarantees atom order)
        molecule = Molecule.from_mapped_smiles("[O:1]([H:2])[H:3]")
        topology = molecule.to_topology().to_openmm()

        # Expected constraints
        constraints = (set(), set(), {(0, 1), (0, 2), (1, 2)})

        for forcefield_path in ("opc3.offxml", "openff_unconstrained-2.3.0.offxml"):
            # Make force field
            forcefield = openmm.app.ForceField()
            forcefield.registerTemplateGenerator(
                SMIRNOFFTemplateGenerator(molecules=molecule, forcefield=forcefield_path).generator
            )

            # We should always get rigid water no matter what is asked for
            assert self.get_terms(forcefield.createSystem(topology, rigidWater=None)) == constraints
            assert self.get_terms(forcefield.createSystem(topology, rigidWater=False)) == constraints
            assert self.get_terms(forcefield.createSystem(topology, rigidWater=True)) == constraints

    def test_constraints_distance(self):
        """
        Tests that constraint distances come from harmonic bonds if created
        through OpenMM and from SMIRNOFF constraints if present.
        """

        from openmm import unit
        from openmm.app import AllBonds

        molecule = Molecule.from_smiles("[F][F]")
        topology = molecule.to_topology().to_openmm()

        offxml = """<?xml version="1.0" encoding="utf-8"?>
<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">{constraints}
    <Bonds version="0.4" potential="harmonic" fractional_bondorder_method="AM1-Wiberg" fractional_bondorder_interpolation="linear">
        <Bond smirks="[#9:1]-[#9:2]" length="2 * angstrom ** 1" k="1000.0 * angstrom ** -2 * kilocalorie_per_mole ** 1"></Bond>
    </Bonds>
    <LibraryCharges version="0.3">
        <LibraryCharge smirks="[#9:1]" charge1="0 * elementary_charge ** 1"></LibraryCharge>
    </LibraryCharges>
    <vdW version="0.4" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom ** 1" switch_width="1.0 * angstrom ** 1" periodic_method="cutoff" nonperiodic_method="no-cutoff">
        <Atom smirks="[#9:1]" epsilon="0 * kilocalorie_per_mole ** 1" sigma="0 * angstrom ** 1"></Atom>
    </vdW>
</SMIRNOFF>"""
        constraints = """
    <Constraints version="0.3">
        <Constraint smirks="[#9:1]-[#9:2]" distance="3 * angstrom ** 1"></Constraint>
    </Constraints>"""

        # Make a system with a SMIRNOFF force field that has only a harmonic
        # bond, and use OpenMM to add a constraint: should get the harmonic
        # bond distance
        generator = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield=offxml.format(constraints=""))
        forcefield = openmm.app.ForceField()
        forcefield.registerTemplateGenerator(generator.generator)
        system = forcefield.createSystem(topology, constraints=AllBonds)
        assert np.isclose(system.getConstraintParameters(0)[2].value_in_unit(unit.angstrom), 2)

        # Now use a SMIRNOFF force field with explicit <Constraints>: even with
        # the harmonic bond present, should get the constraint distance.
        generator = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield=offxml.format(constraints=constraints))
        forcefield = openmm.app.ForceField()
        forcefield.registerTemplateGenerator(generator.generator)
        system = forcefield.createSystem(topology, constraints=AllBonds)
        assert np.isclose(system.getConstraintParameters(0)[2].value_in_unit(unit.angstrom), 3)

    def test_energies_virtual_sites(self):
        """Test potential energies match for systems with virtual sites"""

        test_dir_path = get_data_filename("test_vsites_mols")
        for molecule_path in glob.glob(os.path.join(test_dir_path, "*.sdf")):
            custom_ff = OFFForceField("openff-2.3.0.offxml", get_data_filename("test-virtual-sites.offxml"))

            # These test files have charges that we want to use instead of
            # having these handlers generate them (this is specified also by
            # using charge_from_molecules below, but we delete these to ensure
            # no possibility of other charges being involved)
            custom_ff.deregister_parameter_handler("NAGLCharges")
            custom_ff.deregister_parameter_handler("LibraryCharges")

            # Try non-default values for 1-4 interaction scaling
            custom_ff.get_parameter_handler("vdW").scale14 = 0.7
            custom_ff.get_parameter_handler("Electrostatics").scale14 = 0.3

            # Set up OpenMM ForceField with template generator
            molecules = Molecule.from_file(molecule_path)[:1]
            generator = SMIRNOFFTemplateGenerator(molecules=molecules[0], forcefield=custom_ff.to_string())

            openmm_forcefield = openmm.app.ForceField()
            openmm_forcefield.registerTemplateGenerator(generator.generator)

            # Add virtual sites to OpenMM topology
            openff_topology = Topology.from_molecules(molecules)
            openmm_topology = openff_topology.to_openmm()
            positions = openff_topology.get_positions().to_openmm()
            modeller = openmm.app.Modeller(openmm_topology, positions)
            modeller.addExtraParticles(openmm_forcefield)

            # Make OpenFF-created and ForceField-created systems to compare
            # Note that we've loaded two copies of the same molecule with reversed atom orders.
            # Both have identical partial charges assigned (accounting for the rearrangement)
            # so we pick the first one arbitrarily to be the charge reference.
            smirnoff_system = generator._smirnoff_forcefield.create_openmm_system(
                openff_topology, charge_from_molecules=[molecules[0]]
            )
            openmm_system = openmm_forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

            new_positions = self.propagate_dynamics(modeller.positions, openmm_system)
            self.compare_energies(molecules[0].to_hill_formula(), new_positions, openmm_system, smirnoff_system)
            new_positions = self.propagate_dynamics(new_positions, openmm_system)
            self.compare_energies(molecules[0].to_hill_formula(), new_positions, openmm_system, smirnoff_system)

    def test_energies_multiple_residue(self):
        """Test parameterizing a multi-residue molecule"""

        pdb = PDBFile(get_data_filename("test-ala-3.pdb"))
        molecules = [Molecule.from_topology(Topology.from_pdb(get_data_filename("test-ala-3.pdb")))]
        generator = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield="openff_unconstrained-2.3.0.offxml")
        openmm_forcefield = openmm.app.ForceField()
        openmm_forcefield.registerTemplateGenerator(generator.generator)

        modeller = openmm.app.Modeller(pdb.topology, pdb.positions)
        modeller.addExtraParticles(openmm_forcefield)

        smirnoff_system = generator._smirnoff_forcefield.create_openmm_system(
            Topology.from_openmm(pdb.topology, molecules)
        )
        openmm_system = openmm_forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

        new_positions = self.propagate_dynamics(modeller.positions, openmm_system)
        self.compare_energies("test_energies_multiple_residue", new_positions, openmm_system, smirnoff_system)
        new_positions = self.propagate_dynamics(new_positions, openmm_system)
        self.compare_energies("test_energies_multiple_residue", new_positions, openmm_system, smirnoff_system)

    def test_bespoke_force_field(self):
        """
        Make sure a molecule can be parameterised using a bespoke force field passed as a string to
        the template generator.
        """

        custom_sage = OFFForceField("openff-2.0.0.offxml")
        # Create a simple molecule with one bond type
        ethane = Molecule.from_smiles("CC")
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
        openmm_system = openmm_forcefield.createSystem(
            ethane.to_topology().to_openmm(),
            removeCMMotion=False,
            nonbondedMethod=NoCutoff,
        )

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

        from openmm import unit

        custom_sage = OFFForceField("openff-2.0.0.offxml")
        custom_sage.get_parameter_handler("vdW").scale14 = 0.7
        custom_sage.get_parameter_handler("Electrostatics").scale14 = 0.3

        # Create a simplest 1-4 bond molecule
        ethane = Molecule.from_mapped_smiles("[C:1]([C:2]([H:6])([H:7])[H:8])([H:3])([H:4])[H:5]")

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
            nonbondedMethod=NoCutoff,
        )

        # Find NonbondedForce
        nb_force = [force for force in openmm_system.getForces() if force.__class__.__name__ == "NonbondedForce"][0]

        # Get parameters and exceptions to check
        parameters = []
        for particle_index in range(nb_force.getNumParticles()):
            q, _, epsilon = nb_force.getParticleParameters(particle_index)
            parameters.append(
                (q.value_in_unit(unit.elementary_charge), epsilon.value_in_unit(unit.kilojoule_per_mole))
            )
        exceptions = {}
        for exception_index in range(nb_force.getNumExceptions()):
            i, j, qq, _, epsilon = nb_force.getExceptionParameters(exception_index)
            exceptions[min(i, j), max(i, j)] = (
                qq.value_in_unit(unit.elementary_charge**2),
                epsilon.value_in_unit(unit.kilojoule_per_mole),
            )

        # Expected 1-2 and 1-3 exception pairs (should be zeroed)
        expected_zeroed = {
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (0, 5),
            (0, 6),
            (0, 7),
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (1, 6),
            (1, 7),
            (2, 3),
            (2, 4),
            (3, 4),
            (5, 6),
            (5, 7),
            (6, 7),
        }

        # Expected 1-4 exception pairs (should have requested scales applied)
        expected_scaled = {
            (2, 5),
            (2, 6),
            (2, 7),
            (3, 5),
            (3, 6),
            (3, 7),
            (4, 5),
            (4, 6),
            (4, 7),
        }

        # Check that scaling factors were applied as requested
        assert set(exceptions) == expected_zeroed | expected_scaled
        for i, j in expected_zeroed:
            assert exceptions[i, j] == (0.0, 0.0)
        for i, j in expected_scaled:
            assert np.isclose(exceptions[i, j][0], 0.3 * parameters[i][0] * parameters[j][0])
            assert np.isclose(exceptions[i, j][1], 0.7 * np.sqrt(parameters[i][1] * parameters[j][1]))

    def test_charge_none(self):
        """Test that charges are nonzero after charging if the molecule has None for user charges"""
        for molecule in self.charge_test_molecules:
            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(molecule, None)

            # Molecule should have partial_charges set to None but system should
            # have non-zero charges assigned
            assert molecule.partial_charges is None
            assert not np.allclose(self.charges_from_system(system), 0)

    def test_charge_zero(self):
        """Test that charges are nonzero after charging if the molecule has zero for user charges"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(
                    molecule, unit.Quantity(np.zeros(molecule.n_atoms), unit.elementary_charge)
                )

            # Molecule should have partial_charges set to zero but system should
            # have non-zero charges assigned
            assert molecule.partial_charges is not None
            assert np.all(molecule.partial_charges.m_as(unit.elementary_charge) == 0)
            assert not np.allclose(self.charges_from_system(system), 0)

    def test_charge_valid(self):
        """Test that user-specified partial charges are used if requested"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            user_charges = (
                np.linspace(-1, 1, molecule.n_atoms)
                + molecule.total_charge.m_as(unit.elementary_charge) / molecule.n_atoms
            )

            # Charge mismatch warning should not be raised
            with warnings.catch_warnings():
                warnings.filterwarnings("error", message=CHARGE_WARNING_PATTERN)
                system = self.parameterize_with_charges(molecule, unit.Quantity(user_charges, unit.elementary_charge))

            # User, molecule, and system charges should be equal
            assert molecule.partial_charges is not None
            assert self.charges_are_equal(system, molecule)
            assert np.allclose(self.charges_from_system(system), user_charges)

    def test_charge_warning(self):
        """Test that a warning is raised with user-specified partial charges that do not sum to the total charge"""
        from openff.units import unit

        for molecule in self.charge_test_molecules:
            user_charges = (
                np.linspace(-1, 1, molecule.n_atoms)
                + molecule.total_charge.m_as(unit.elementary_charge) / molecule.n_atoms
                + 1
            )

            # Charge mismatch warning should be raised
            with pytest.warns(match=CHARGE_WARNING_PATTERN):
                system = self.parameterize_with_charges(molecule, unit.Quantity(user_charges, unit.elementary_charge))

            # User, molecule, and system charges should be equal
            assert molecule.partial_charges is not None
            assert self.charges_are_equal(system, molecule)
            assert np.allclose(self.charges_from_system(system), user_charges)


@pytest.mark.espaloma
class TestEspalomaTemplateGenerator(TemplateGeneratorBaseCase):
    TEMPLATE_GENERATOR = EspalomaTemplateGenerator

    def test_retrieve_forcefields(self):
        """Test a force field can be retrieved"""
        # Test loading model by specifying version number
        generator = EspalomaTemplateGenerator(forcefield="espaloma-0.3.2")
        del generator
        # Test loading model from remote URL
        url = "https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt"
        generator = EspalomaTemplateGenerator(forcefield=url)
        del generator
        # Test loading model from filename
        with tempfile.TemporaryDirectory() as tmpdirname:
            import urllib

            filename = os.path.join(tmpdirname, "model.pt")
            urllib.request.urlretrieve(url, filename=filename)
            # Create a new database file
            generator = EspalomaTemplateGenerator(forcefield=filename)
            del generator

    def test_energies(self):
        """Test potential energies match between openff-toolkit and OpenMM ForceField"""

        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in EspalomaTemplateGenerator.INSTALLED_FORCEFIELDS:
            print(f"Testing energies for {small_molecule_forcefield}...")
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
                openmm_system = openmm_forcefield.createSystem(
                    molecule.to_topology().to_openmm(),
                    removeCMMotion=False,
                    nonbondedMethod=NoCutoff,
                )

                # Retrieve System generated by Espaloma
                espaloma_system = generator.get_openmm_system(molecule)

                # Compare energies and forces
                self.compare_energies_single(molecule, openmm_system, espaloma_system)

                # Run some dynamics
                molecule = self.propagate_dynamics_single(molecule, espaloma_system)

                # Compare energies again
                self.compare_energies_single(molecule, openmm_system, espaloma_system)

    def test_partial_charges_are_none(self):
        """Test parameterizing a small molecule with `partial_charges=None` instead
        of zeros (happens frequently in OFFTK>=0.7.0)"""
        from openff.units import unit

        molecule = Molecule.from_smiles("C=O")
        molecule.generate_conformers(n_conformers=1)
        # molecule._partial_charges = None
        assert (molecule.partial_charges is None) or np.all(molecule.partial_charges / unit.elementary_charge == 0)
        # Test all supported SMIRNOFF force fields
        for small_molecule_forcefield in EspalomaTemplateGenerator.INSTALLED_FORCEFIELDS:
            print(f"Testing energies for {small_molecule_forcefield}...")
            # Create a generator that knows about a few molecules
            # TODO: Should the generator also load the appropriate force field files into the ForceField object?
            generator = EspalomaTemplateGenerator(molecules=[molecule], forcefield=small_molecule_forcefield)
            # Create a ForceField
            openmm_forcefield = openmm.app.ForceField()
            # Register the template generator
            openmm_forcefield.registerTemplateGenerator(generator.generator)
            # Create OpenMM System using OpenMM app
            openmm_forcefield.createSystem(
                molecule.to_topology().to_openmm(),
                removeCMMotion=False,
                nonbondedMethod=NoCutoff,
            )
            generator.get_openmm_system(molecule)

    # def test_template_generator_kwargs(self):
    #     """Test """

    def test_keyword_arguments_default(self):
        """
        Test the default behavior for the keyword arguments, which is using "openff_unconstrained-2.2.1"
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
        # Assert the reference forcefield is the default "openff_unconstrained-2.2.1"
        default_ref_ff = "openff_unconstrained-2.2.1"
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
        generator = EspalomaTemplateGenerator(
            molecules=[molecule],
            forcefield="espaloma-0.3.2",
            template_generator_kwargs=espaloma_generator_kwargs,
        )
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
