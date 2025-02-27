import argparse
import copy
import enum
import itertools
import numpy
import openmm.app
import os
import subprocess
import sys
import tempfile
import yaml

# The CHARMM value can be found throughout the CHARMM source code.  The OpenMM
# value is from SimTKOpenMMRealType.h.
CHARMM_VACUUM_PERMITTIVITY = (
    openmm.unit.AVOGADRO_CONSTANT_NA
    * openmm.unit.elementary_charge**2
    / (4 * numpy.pi * 332.0716 * openmm.unit.kilocalorie_per_mole * openmm.unit.angstrom)
)
OPENMM_VACUUM_PERMITTIVITY = 8.8541878128e-12 * openmm.unit.farad / openmm.unit.meter
CHARMM_ELECTROSTATIC_SCALE = CHARMM_VACUUM_PERMITTIVITY / OPENMM_VACUUM_PERMITTIVITY

DEFAULT_ABSOLUTE_ENERGY_TOLERANCE = 2e-4  # kcal/mol
DEFAULT_RELATIVE_ENERGY_TOLERANCE = 2e-4
DEFAULT_ABSOLUTE_FORCE_TOLERANCE = 2e-2  # kcal/mol/Å
DEFAULT_RELATIVE_FORCE_TOLERANCE = 2e-2

DEFAULT_PERTURB_DISTANCE = 0.05  # Å

DEFAULT_OPENMM_PLATFORM = "automatic"

DEFAULT_TOPPAR_DIRECTORY = "tests/toppar"
DEFAULT_FFXML_DIRECTORY = "ffxml"

ENERGY_PRINT_UNIT = openmm.unit.kilocalorie_per_mole
FORCE_PRINT_UNIT = ENERGY_PRINT_UNIT / openmm.unit.angstrom

# A string that shouldn't appear in the CHARMM output.
ENERGY_FORCE_TAG = "OPENMMFORCEFIELDS-TEST-CHARMM-ENERGY-FORCES"

CHARMM_ENERGY_CHECK_PRECISION = 1e-3

OPENMM_CONSTRAINTS = {
    "None": None,
    "HBonds": openmm.app.HBonds,
    "AllBonds": openmm.app.AllBonds,
    "HAngles": openmm.app.HAngles,
}

OPENMM_NONBONDED_METHOD = {
    "NoCutoff": openmm.app.NoCutoff,
    "CutoffNonPeriodic": openmm.app.CutoffNonPeriodic,
    "CutoffPeriodic": openmm.app.CutoffPeriodic,
    "Ewald": openmm.app.Ewald,
    "PME": openmm.app.PME,
    "LJPME": openmm.app.LJPME,
}

DO_COLOR = sys.stdout.isatty()


class ForceGroup(enum.Enum):
    """
    Force groups for decomposing contributions to the total energy.  Bonds and
    Urey-Bradley forces are combined since these can both be represented in
    OpenMM in a single HarmonicBondForce, making them impossible to isolate from
    each other.
    """

    BONDS_UREY_BRADLEY = 1
    ANGLES = 2
    PROPERS = 3
    IMPROPERS = 4
    CMAP = 5
    NONBONDED = 6
    DRUDE = 7


# For each force group, the CHARMM skip keywords, the labels in the CHARMM
# output, and the OpenMM force types.
FORCE_GROUP_DATA = {
    ForceGroup.BONDS_UREY_BRADLEY: (("bond", "urey"), ("BONDs", "UREY-b"), (openmm.HarmonicBondForce,)),
    ForceGroup.ANGLES: (("angl",), ("ANGLes",), (openmm.HarmonicAngleForce,)),
    ForceGroup.PROPERS: (("dihe",), ("DIHEdrals",), (openmm.PeriodicTorsionForce,)),
    ForceGroup.IMPROPERS: (("impr",), ("IMPRopers",), (openmm.CustomTorsionForce,)),
    ForceGroup.CMAP: (("cmap",), ("CMAPs",), (openmm.CMAPTorsionForce,)),
    ForceGroup.NONBONDED: (
        ("vdw", "elec", "imnb", "imel", "ewks", "ewse", "ewex"),
        ("VDWaals", "ELEC", "IMNBvdw", "IMELec", "EWKSum", "EWSElf", "EWEXcl"),
        (openmm.NonbondedForce, openmm.CustomNonbondedForce, openmm.CustomBondForce),
    ),
    ForceGroup.DRUDE: (
        (),
        (),
        (openmm.DrudeForce,),
    ),
    # NOTE: Impropers could be under a PeriodicTorsionForce but this almost
    # never occurs in the CHARMM force field.  The decomposition between propers
    # and impropers would fail in that case even if the total energy was being
    # evaluated correctly, but this shouldn't happen for most systems.
}

# Groups to use when evaluating just electrostatic forces (to scale CHARMM vs.
# OpenMM results accurately).
CHARMM_ELECTROSTATIC_GROUPS = ("elec", "imel", "ewks", "ewse", "ewex")

# Groups to merge into Drude forces when comparing CHARMM with OpenMM (since
# CHARMM uses a different splitting of Drude forces).
CHARMM_DRUDE_GROUPS = (
    ForceGroup.BONDS_UREY_BRADLEY,
    ForceGroup.NONBONDED,
)


def main():
    """
    Main entry point for test script.  Parses command-line arguments and runs.
    """

    parser = argparse.ArgumentParser(
        description="Test script for CHARMM force field conversion to OpenMM FFXML format"
    )

    parser.add_argument(
        "test_path",
        nargs="*",
        help="Test specification (in a YAML file) to run",
    )

    # Flags to compute energies and forces using particular methods.
    parser.add_argument(
        "--charmm",
        action="store_true",
        help="Compute energies and forces with CHARMM (requires charmm executable in PATH)",
    )
    parser.add_argument(
        "--openmm-charmm",
        action="store_true",
        help="Compute energies and forces with OpenMM using OpenMM CharmmParameterSet",
    )
    parser.add_argument(
        "--parmed-charmm",
        action="store_true",
        help="Compute energies and forces with OpenMM using ParmEd CharmmParameterSet",
    )
    parser.add_argument(
        "--openmm-ffxml",
        action="store_true",
        help="Compute energies and forces with OpenMM using OpenMM ForceField",
    )

    # Flags to control energy and force comparison.
    parser.add_argument(
        "--skip-term",
        action="append",
        default=[],
        choices=[force_group.name for force_group in ForceGroup],
        help="Exclude the selected force group from the energies considered",
    )
    parser.add_argument(
        "--absolute-energy-tolerance",
        type=float,
        default=DEFAULT_ABSOLUTE_ENERGY_TOLERANCE,
        metavar="kcal/mol",
        help=f"Absolute tolerance for matching energies (default {DEFAULT_ABSOLUTE_ENERGY_TOLERANCE} kcal/mol)",
    )
    parser.add_argument(
        "--relative-energy-tolerance",
        type=float,
        default=DEFAULT_RELATIVE_ENERGY_TOLERANCE,
        metavar="ratio",
        help=f"Relative tolerance for matching energies (default {DEFAULT_RELATIVE_ENERGY_TOLERANCE})",
    )
    parser.add_argument(
        "--absolute-force-tolerance",
        type=float,
        default=DEFAULT_ABSOLUTE_FORCE_TOLERANCE,
        metavar="kcal/mol/Å",
        help=f"Absolute tolerance for matching forces (default {DEFAULT_ABSOLUTE_FORCE_TOLERANCE} kcal/mol/Å)",
    )
    parser.add_argument(
        "--relative-force-tolerance",
        type=float,
        default=DEFAULT_RELATIVE_FORCE_TOLERANCE,
        metavar="ratio",
        help=f"Relative tolerance for matching forces (default {DEFAULT_RELATIVE_FORCE_TOLERANCE})",
    )

    # Flags to control adding replicates involving random coordinate
    # perturbation.
    parser.add_argument(
        "--perturb-replicates",
        type=int,
        default=0,
        metavar="count",
        help="Repeat tests with perturbed positions this many times",
    )
    parser.add_argument(
        "--perturb-distance",
        type=float,
        default=DEFAULT_PERTURB_DISTANCE,
        metavar="Å",
        help=f"Perturb positions by up to this distance (default {DEFAULT_PERTURB_DISTANCE} Å)",
    )
    parser.add_argument(
        "--perturb-seed",
        type=int,
        metavar="seed",
        help="Random seed for repeating tests with perturbed positions",
    )

    # Miscellaneous options.
    parser.add_argument(
        "--dump",
        action="store_true",
        help="Dump input files for debugging",
    )
    parser.add_argument(
        "--openmm-platform",
        default=DEFAULT_OPENMM_PLATFORM,
        choices=["automatic"]
        + [openmm.Platform.getPlatform(index).getName() for index in range(openmm.Platform.getNumPlatforms())],
        help='Platform for OpenMM to use, or "automatic" to let OpenMM choose (default {DEFAULT_OPENMM_PLATFORM})',
    )
    parser.add_argument(
        "--openmm-ffxml-fix-impropers",
        action="store_true",
        help="Try to repair the ordering of impropers in OpenMM FFXML-created systems",
    )

    # Specify directories where parameter files can be found.
    parser.add_argument(
        "--toppar-directory",
        default=DEFAULT_TOPPAR_DIRECTORY,
        metavar="path",
        help="Path to CHARMM parameter files (default toppar)",
    )
    parser.add_argument(
        "--ffxml-directory",
        default=DEFAULT_FFXML_DIRECTORY,
        metavar="path",
        help="Path to OpenMM force field files (default ffxml)",
    )

    run(parser.parse_args())


def run(args):
    """
    Runs the specified tests with the specified options and exits.  The return
    code will be the number of failed tests.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line arguments.
    """

    if args.perturb_replicates:
        # Passing None (when --perturb-seed is not given) indicates to NumPy
        # that a random seed should be generated automatically.
        perturb_spec = (
            numpy.random.default_rng(seed=args.perturb_seed),
            args.perturb_replicates,
            args.perturb_distance,
        )
    else:
        perturb_spec = None

    skip_name_set = set(args.skip_term)
    skip_set = set(force_group for force_group in ForceGroup if force_group.name in skip_name_set)

    runner = TestRunner(
        args.charmm,
        args.openmm_charmm,
        args.parmed_charmm,
        args.openmm_ffxml,
        skip_set,
        args.absolute_energy_tolerance * openmm.unit.kilocalorie_per_mole,
        args.relative_energy_tolerance,
        args.absolute_force_tolerance * openmm.unit.kilocalorie_per_mole / openmm.unit.angstrom,
        args.relative_force_tolerance,
        perturb_spec,
        args.dump,
        args.openmm_platform,
        args.openmm_ffxml_fix_impropers,
        args.toppar_directory,
        args.ffxml_directory,
    )

    test_results = []
    for test_path in args.test_path:
        try:
            term_failure_count = runner.run_test(test_path)
            test_results.append((test_path, term_failure_count == 0, f"{term_failure_count} terms exceeded tolerance"))
        except Exception as error:
            test_results.append((test_path, False, f"{type(error).__name__}: {error}"))

    if not test_results:
        print("No tests run (use --help for information on how to specify tests)")
        sys.exit()

    # Summarize test results.
    test_count = len(test_results)
    test_names, success_flags, _ = zip(*test_results)
    success_count = sum(success_flags)
    failure_count = test_count - success_count

    print(f"Ran {test_count} tests ({success_count} successes, {failure_count} failures)")
    test_name_width = max(map(len, test_names))
    for test_name, success_flag, test_message in test_results:
        formatted_line = f"{'Succeeded' if success_flag else 'Failed':9}  {test_message}"
        print(
            f"    {test_name:{test_name_width}}  {formatted_line if success_flag else color_message(formatted_line)}"
        )

    sys.exit(failure_count)


class TestRunner:
    """
    Holds options used for running tests.

    Parameters
    ----------
    do_charmm : bool
        Whether or not to compute energies and forces with CHARMM.
    do_openmm_charmm : bool
        Whether or not to compute energies and forces with OpenMM using OpenMM
        CharmmParameterSet.
    do_parmed_charmm : bool
        Whether or not to compute energies and forces with OpenMM using ParmEd
        CharmmParameterSet.
    do_openmm_ffxml : bool
        Whether or not to compute energies and forces with OpenMM using OpenMM
        ForceField.
    skip_set : set(ForceGroup)
        Set of terms to ignore during testing.
    absolute_energy_tolerance : openmm.unit.Quantity
        Absolute tolerance for matching energies.
    relative_energy_tolerance : float
        Relative tolerance for matching energies.
    absolute_force_tolerance : openmm.unit.Quantity
        Absolute tolerance for matching forces.
    relative_force_tolerance : float
        Relative tolerance for matching forces.
    perturb_spec : (numpy.random.Generator, int, float) or None
        If not None, repeat tests the given number of times with positions
        perturbed by up to the given distance in angstroms using random numbers
        from the given generator.
    do_dump : bool
        Whether or not to dump input files for debugging.
    openmm_platform : str
        Name of an OpenMM Platform to use, or "automatic" to let OpenMM choose.
    openmm_ffxml_fix_impropers : bool
        Whether or not to try to repair the ordering of impropers in OpenMM
        FFXML-created systems.
    toppar_directory : str
        Directory containing CHARMM parameter files.
    ffxml_directory : str
        Directory containing OpenMM force field files.
    """

    def __init__(
        self,
        do_charmm,
        do_openmm_charmm,
        do_parmed_charmm,
        do_openmm_ffxml,
        skip_set,
        absolute_energy_tolerance,
        relative_energy_tolerance,
        absolute_force_tolerance,
        relative_force_tolerance,
        perturb_spec,
        do_dump,
        openmm_platform,
        openmm_ffxml_fix_impropers,
        toppar_directory,
        ffxml_directory,
    ):
        self.do_charmm = do_charmm
        self.do_openmm_charmm = do_openmm_charmm
        self.do_parmed_charmm = do_parmed_charmm
        self.do_openmm_ffxml = do_openmm_ffxml
        self.skip_set = skip_set
        self.absolute_energy_tolerance = absolute_energy_tolerance
        self.relative_energy_tolerance = relative_energy_tolerance
        self.absolute_force_tolerance = absolute_force_tolerance
        self.relative_force_tolerance = relative_force_tolerance
        self.perturb_spec = perturb_spec
        self.do_dump = do_dump
        self.openmm_platform = openmm_platform
        self.openmm_ffxml_fix_impropers = openmm_ffxml_fix_impropers
        self.toppar_directory = toppar_directory
        self.ffxml_directory = ffxml_directory

    def run_test(self, test_path):
        """
        Runs the specified test with the stored options.

        Parameters
        ----------
        test_path : str
            Path to a YAML file containing a test specification.

        Returns
        -------
        int
            The number of failed checks.
        """

        test_directory = os.path.dirname(test_path)
        with open(test_path) as test_spec_file:
            test_spec = yaml.safe_load(test_spec_file)

        # Insert default values into test specification.
        for default_list in ("charmm_files", "ffxml_files", "charmm_commands"):
            if default_list not in test_spec:
                test_spec[default_list] = []

        # Adjust system options read from specification.
        if "box" in test_spec:
            box = test_spec["box"]
            for key in ("a", "b", "c"):
                box[key] = yaml_to_quantity(box[key])
        if "create_system_options" in test_spec:
            create_system_options = test_spec["create_system_options"]
            for key, value in create_system_options.items():
                if key in ("flexibleConstraints", "rigidWater", "removeCMMotion", "useDispersionCorrection"):
                    # Boolean values can be read directly from the test file.
                    pass
                elif key in ("drudeMass", "hydrogenMass", "nonbondedCutoff", "switchDistance"):
                    create_system_options[key] = yaml_to_quantity(value)
                elif key == "constraints":
                    create_system_options[key] = OPENMM_CONSTRAINTS[value]
                elif key == "nonbondedMethod":
                    create_system_options[key] = OPENMM_NONBONDED_METHOD[value]
                else:
                    raise ValueError(f"unsupported system option {key!r}")
        else:
            test_spec["create_system_options"] = {}

        # Generate reference coordinates.
        coordinates_list = [openmm.app.PDBFile(os.path.join(test_directory, test_spec["pdb_file"])).getPositions()]
        if self.perturb_spec is not None:
            generator, perturb_replicates, perturb_distance = self.perturb_spec
            for replicate in range(perturb_replicates):
                perturbed_coordinates = coordinates_list[0].value_in_unit(openmm.unit.angstrom).copy()
                for index in range(len(perturbed_coordinates)):
                    perturbed_coordinates[index] += openmm.Vec3(*(perturb_distance * random_in_unit_ball(generator)))
                coordinates_list.append(openmm.unit.Quantity(perturbed_coordinates, openmm.unit.angstrom))

        dump_path = test_path.replace("/", "_")

        # Run tests and get results.
        results_sets = {}
        if self.do_charmm:
            results_sets["charmm"] = self.get_charmm_energies_forces(
                test_directory,
                test_spec,
                coordinates_list,
                f"{dump_path}_charmm.inp" if self.do_dump else None,
                f"{dump_path}_charmm.log" if self.do_dump else None,
            )
        if self.do_openmm_charmm:
            results_sets["openmm_charmm"] = self.get_openmm_charmm_energies_forces(
                test_directory,
                test_spec,
                coordinates_list,
                f"{dump_path}_openmm_charmm.xml" if self.do_dump else None,
            )
        if self.do_parmed_charmm:
            results_sets["parmed_charmm"] = self.get_parmed_charmm_energies_forces(
                test_directory,
                test_spec,
                coordinates_list,
                f"{dump_path}_parmed_charmm.xml" if self.do_dump else None,
            )
        if self.do_openmm_ffxml:
            results_sets["openmm_ffxml"] = self.get_openmm_ffxml_energies_forces(
                test_directory,
                test_spec,
                coordinates_list,
                f"{dump_path}_openmm_ffxml.xml" if self.do_dump else None,
            )

        # Zero anything that's not being tested.
        for results in results_sets.values():
            for result in results:
                for skip_group in self.skip_set:
                    result[None]["energy"] -= result[skip_group]["energy"]
                    result[skip_group]["energy"] *= 0
                    result[None]["forces"] -= result[skip_group]["forces"]
                    result[skip_group]["forces"] *= 0

        if len(results_sets) < 2:
            raise ValueError("must specify at least 2 methods to compare")

        # Special handling occurs for Drude systems.
        is_drude = test_spec.get("drude", False)

        print(f"Test {test_spec['name']!r} in {test_path!r}:")
        failure_count = 0
        for (name_1, results_1), (name_2, results_2) in itertools.combinations(results_sets.items(), 2):
            print(f"    Comparing {name_1} vs. {name_2}")

            # If comparing Drude energies and forces between CHARMM and OpenMM,
            # lump bonds, non-bonded forces, and Drude forces together (here
            # into the DRUDE group) since they are split up differently between
            # the two programs.
            if is_drude and name_1 == "charmm" and name_2 != "charmm":
                results_1 = copy.deepcopy(results_1)
                results_2 = copy.deepcopy(results_2)
                for results in (results_1, results_2):
                    for result in results:
                        for group in CHARMM_DRUDE_GROUPS:
                            result[ForceGroup.DRUDE]["energy"] += result[group]["energy"]
                            result[group]["energy"] *= 0
                            result[ForceGroup.DRUDE]["forces"] += result[group]["forces"]
                            result[group]["forces"] *= 0

            # Print results for energies.
            print(f"        {'Energy error':20}  {'|ΔE| (kcal/mol)':20}  {'Relative':12}")
            for force_group in (None, *ForceGroup):
                difference_data = []
                for result_1, result_2 in zip(results_1, results_2):
                    energy_1 = result_1[force_group]["energy"]
                    energy_2 = result_2[force_group]["energy"]
                    difference_data.append((energy_1, energy_2, *self._compare_energies(energy_1, energy_2)))
                failure_count += self._format_energy_difference_data(
                    "TOTAL" if force_group is None else force_group.name, difference_data
                )
            print()

            # Print results for forces.
            print(f"        {'Force error (max)':20}  {'|Δf| (kcal/mol/Å)':20}  {'Relative':12}")
            for force_group in (None, *ForceGroup):
                difference_data = []
                for result_1, result_2 in zip(results_1, results_2):
                    forces_1 = result_1[force_group]["forces"]
                    forces_2 = result_2[force_group]["forces"]
                    is_virtual_1 = result_1[force_group]["mask"] != 0
                    is_virtual_2 = result_2[force_group]["mask"] != 0

                    if not numpy.all(is_virtual_1 == is_virtual_2):
                        raise ValueError("inconsistent specification of virtual sites")

                    forces_1[is_virtual_1] *= 0
                    forces_2[is_virtual_2] *= 0

                    magnitudes_1 = numpy.linalg.norm(forces_1.value_in_unit(FORCE_PRINT_UNIT), axis=1)
                    magnitudes_2 = numpy.linalg.norm(forces_2.value_in_unit(FORCE_PRINT_UNIT), axis=1)
                    absolute_differences = numpy.linalg.norm(
                        (forces_2 - forces_1).value_in_unit(FORCE_PRINT_UNIT), axis=1
                    )
                    relative_mask = (magnitudes_1 != 0) | (magnitudes_2 != 0)
                    relative_differences = numpy.zeros(absolute_differences.shape)
                    relative_differences[relative_mask] = absolute_differences[relative_mask] / numpy.maximum(
                        magnitudes_1[relative_mask], magnitudes_2[relative_mask]
                    )
                    absolute_failures = absolute_differences > self.absolute_force_tolerance.value_in_unit(
                        FORCE_PRINT_UNIT
                    )
                    relative_failures = relative_differences > self.relative_force_tolerance
                    absolute_differences *= FORCE_PRINT_UNIT

                    difference_data.append(
                        (absolute_failures, relative_failures, absolute_differences, relative_differences)
                    )

                failure_count += self._format_force_difference_data(
                    "TOTAL" if force_group is None else force_group.name, difference_data
                )
            print()

        return failure_count

    def _compare_energies(self, energy_1, energy_2):
        """
        Helper function for comparing energies.

        Parameters
        ----------
        energy_1 : openmm.unit.Quantity
            The first energy to compare.
        energy_2 : openmm.unit.Quantity
            The second energy to compare.

        Returns
        -------
        openmm.unit.Quantity, float
            The magnitude of the difference between energies, and the relative
            difference.  Both values will be non-negative.  The relative
            difference will use whatever energy value has the larger magnitude
            as a reference scale, and will be zero if both energies are zero.
        """

        absolute_difference = abs(energy_2 - energy_1)
        if energy_1 or energy_2:
            relative_difference = absolute_difference / max(abs(energy_1), abs(energy_2))
        else:
            relative_difference = 0.0
        return absolute_difference, relative_difference

    def _format_energy_difference_data(self, label, difference_data):
        """
        Prints a table of energy difference data.

        Parameters
        ----------
        label : str
            A label for the relevant force group.
        difference_data : list(tuple(openmm.unit.Quantity, openmm.unit.Quantity, openmm.unit.Quantity, float))
            For each set of coordinates tested, energies evaluated in two ways,
            followed by the output of _compare_energies().

        Returns
        -------
        bool
            Whether or not any test failures were detected.
        """

        energies_1, energies_2, absolute_differences, relative_differences = zip(*difference_data)
        failures = [
            absolute_difference > self.absolute_energy_tolerance
            and relative_difference > self.relative_energy_tolerance
            for absolute_difference, relative_difference in zip(absolute_differences, relative_differences)
        ]

        any_failures = any(failures)
        formatted_line = (
            f"{max(absolute_differences).value_in_unit(ENERGY_PRINT_UNIT):20.12e}  {max(relative_differences):12.8f}"
        )
        if any_failures:
            formatted_line = f"{color_message(formatted_line)}  E_1 (kcal/mol)        E_2 (kcal/mol)"
        print(f"        {label:20}  {formatted_line}")

        if any_failures:
            for energy_1, energy_2, absolute_difference, relative_difference, failure in zip(
                energies_1, energies_2, absolute_differences, relative_differences, failures
            ):
                formatted_line = (
                    f"{energy_1.value_in_unit(ENERGY_PRINT_UNIT):20.12e}  "
                    f"{energy_2.value_in_unit(ENERGY_PRINT_UNIT):20.12e}"
                )
                print(f"{' ' * 66}{color_message(formatted_line) if failure else formatted_line}")

        return any_failures

    def _format_force_difference_data(self, label, difference_data):
        """
        Prints a table of force difference data.

        Parameters
        ----------
        label : str
            A label for the relevant force group.
        difference_data : list(tuple(bool, bool, openmm.unit.Quantity, float))
            For each set of coordinates tested, for each particle, the output of
            _compare_forces().
        """

        absolute_failures, relative_failures, absolute_differences, relative_differences = zip(*difference_data)

        any_failures = any(
            numpy.any(absolute_failure & relative_failure)
            for absolute_failure, relative_failure in zip(absolute_failures, relative_failures)
        )
        absolute_display = max(
            numpy.amax(absolute_difference.value_in_unit(FORCE_PRINT_UNIT))
            for absolute_difference in absolute_differences
        )
        relative_display = numpy.amax(relative_differences)
        formatted_line = f"{absolute_display:20.12e}  {relative_display:12.8f}"
        if any_failures:
            formatted_line = color_message(formatted_line)
        print(f"        {label:20}  {formatted_line}")

        return any_failures

    def get_charmm_energies_forces(self, test_directory, test_spec, coordinates_list, dump_in_path, dump_out_path):
        """
        Computes energies and forces with CHARMM.

        Parameters
        ----------
        test_directory : str
            Path to the directory containing a test specification YAML file.
        test_spec : dict
            A test specification read from a YAML file.
        coordinates_list : list
            A list of arrays of positions.  Energies and forces will be evaluated
            for each set of coordinates.
        dump_in_path : str or None
            Path to dump input files for debugging, or None to skip dumping.
        dump_out_path : str or None
            Path to dump output files for debugging, or None to skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary of dictionaries for each
            ForceGroup (or None for all force field terms) containing the system
            energy as "energy" and forces as "forces".
        """

        # Prepare the CHARMM input line by line.
        charmm_input_lines = ["* Test", "*"]

        # Load all input files.
        first_rtf = True
        first_para = True
        for charmm_file in test_spec["charmm_files"]:
            if charmm_file.endswith(".rtf") or charmm_file.endswith(".top"):
                append = "" if first_rtf else "append "
                if first_rtf:
                    first_rtf = False
                charmm_input_lines.append(
                    f"read rtf card {append}name {os.path.join(self.toppar_directory, charmm_file)}"
                )
            elif charmm_file.endswith(".prm") or charmm_file.endswith(".par"):
                append = "" if first_para else "append "
                if first_para:
                    first_para = False
                charmm_input_lines.append(
                    f"read para card {append}name {os.path.join(self.toppar_directory, charmm_file)} flex"
                )
            else:
                charmm_input_lines.append(f"stream {os.path.join(self.toppar_directory, charmm_file)}")
        charmm_input_lines.append(f"read psf card name {os.path.join(test_directory, test_spec['psf_file'])}")

        # Execute custom commands.
        charmm_input_lines.extend(test_spec["charmm_commands"])

        with tempfile.TemporaryDirectory() as temporary_path:
            # The PSF file is only used to check for the presence of lone pairs
            # and to retrieve data needed to write the coordinate files.
            psf_file = openmm.app.CharmmPsfFile(os.path.join(test_directory, test_spec["psf_file"]))
            virtual_site_mask = numpy.zeros(len(psf_file.atom_list), dtype=bool)
            for lone_pair in psf_file.lonepair_list:
                virtual_site_mask[lone_pair[0]] = True

            for coordinates_index, coordinates in enumerate(coordinates_list):
                # Load coordinates.  Ignore residue and atom name checks and
                # read sequentially.
                coordinates_path = os.path.join(temporary_path, f"{coordinates_index}.crd")
                charmm_input_lines.append(f"read coor ignore name {coordinates_path}")

                # Set up periodic images.
                if not coordinates_index and "box" in test_spec:
                    box = test_spec["box"]
                    box_a = box["a"].value_in_unit(openmm.unit.angstrom)
                    box_b = box["b"].value_in_unit(openmm.unit.angstrom)
                    box_c = box["c"].value_in_unit(openmm.unit.angstrom)
                    charmm_input_lines.append(f"crys define orth {box_a} {box_b} {box_c} 90 90 90")
                    charmm_input_lines.append(f"crys build cutoff {max(box_a, box_b, box_c)} noper 0")

                # Evaluate the energy and forces for the entire system.  Append
                # a tag that we can search for in the output.
                charmm_input_lines.append(f"!{ENERGY_FORCE_TAG}")
                charmm_input_lines.append("skip none")
                charmm_input_lines.append("ener")
                charmm_input_lines.append("coor force comp")
                charmm_input_lines.append("print coor comp")

                # Evaluate energies and forces for each force group.  Turn off
                # all terms (skip all), but leave (excl) some in.
                for force_group in ForceGroup:
                    charmm_input_lines.append(f"! {ENERGY_FORCE_TAG}")
                    charmm_input_lines.append(f"skip all excl {' '.join(FORCE_GROUP_DATA[force_group][0])}")
                    charmm_input_lines.append("ener")
                    charmm_input_lines.append("coor force comp")
                    charmm_input_lines.append("print coor comp")

                # Evaluate just electrostatic energies and forces.  This is used
                # to get a good comparison between OpenMM and CHARMM because of
                # their different values for the vacuum permittivity.
                charmm_input_lines.append(f"! {ENERGY_FORCE_TAG}")
                charmm_input_lines.append(f"skip all excl {' '.join(CHARMM_ELECTROSTATIC_GROUPS)}")
                charmm_input_lines.append("ener")
                charmm_input_lines.append("coor force comp")
                charmm_input_lines.append("print coor comp")

            charmm_input_lines.append("stop")

            if dump_in_path is not None:
                with open(dump_in_path, "w") as dump_in_file:
                    for line in charmm_input_lines:
                        print(line, file=dump_in_file)

            # Write coordinates into the temporary directory.  Only write X, Y,
            # and Z coordinates and other numeric fields as names will be
            # ignored.  Note that this is a fixed-width format.
            for coordinates_index, coordinates in enumerate(coordinates_list):
                with open(os.path.join(temporary_path, f"{coordinates_index}.crd"), "w") as coordinates_file:
                    coordinates_file.write(f"* TEST\n*\n{len(coordinates):10}  EXT\n")
                    for atom_index, (atom, xyz) in enumerate(zip(psf_file.atom_list, coordinates)):
                        x, y, z = xyz.value_in_unit(openmm.unit.angstrom)
                        coordinates_file.write(
                            f"{atom_index + 1:10}{atom.residue.idx:10}                    "
                            f"{x:20.10f}{y:20.10f}{z:20.10f}            {atom.residue.idx:8}{1:20.10f}\n"
                        )

            # Run CHARMM.
            result = subprocess.run(["charmm"], input="\n".join(charmm_input_lines).encode(), capture_output=True)

        if dump_out_path is not None:
            with open(dump_out_path, "wb") as dump_out_file:
                dump_out_file.write(result.stdout)
                dump_out_file.write(result.stderr)

        # If something went wrong, print the output for debugging.
        if result.returncode:
            print(result.stdout.decode())
            print(result.stderr.decode())
            raise RuntimeError(f"CHARMM exited with code {result.returncode}")

        # Parse CHARMM output.
        lines = iter(result.stdout.decode().split("\n"))

        # Helper function for skipping lines until one matches a predicate.
        def advance_until(predicate):
            while True:
                next_line = next(lines)
                if predicate(next_line):
                    return next_line
            raise ValueError("can't find expected pattern in CHARMM output")

        # Helper function for finding lines until one matches a predicate.
        def yield_until(predicate):
            while True:
                next_line = next(lines)
                if predicate(next_line):
                    return
                yield next_line

        # Helper function for reading energies and forces.
        def parse_energies_forces(force_group):
            # Look for the tag we embedded in the input.
            advance_until(lambda line: ENERGY_FORCE_TAG in line)

            # Extract the energy labels and the values themselves.
            ener_labels = [advance_until(lambda line: line.startswith("ENER "))]
            ener_labels.extend(yield_until(lambda line: "-" * 9 in line))
            ener_values = list(yield_until(lambda line: "-" * 9 in line))

            ener_labels = [labels.strip().split() for labels in ener_labels]

            def split_values(values):
                yield values[:14]
                values = values[14:]
                while values:
                    yield values[:13]
                    values = values[13:]

            ener_values = [list(split_values(values)) for values in ener_values]

            # Helper function for finding the indices associated with a label.
            def find_label(desired_label):
                for index_1, labels in enumerate(ener_labels):
                    for index_2, label in enumerate(labels):
                        if label == desired_label:
                            return index_1, index_2
                return None

            # Helper function for finding an energy value.
            def find_value(desired_label):
                # If there is no force (e.g., CMAP for a monomer) in the system,
                # CHARMM may not write it out.
                indices = find_label(desired_label)
                if indices is None:
                    return 0.0
                index_1, index_2 = indices

                # Indices must be adjusted for correct parsing.
                index_2 -= 1
                if not index_1:
                    index_2 -= 1
                return float(ener_values[index_1][index_2])

            # Get the total energy.
            energy = find_value("ENERgy")

            # Check to make sure that the energy decomposition makes sense.
            expected_categories = set(
                category
                for check_force_group in (ForceGroup if force_group is None else (force_group,))
                for category in FORCE_GROUP_DATA[check_force_group][1]
            )
            expected_energy = 0.0
            for check_force_group in ForceGroup:
                for category in FORCE_GROUP_DATA[check_force_group][1]:
                    category_energy = find_value(category)
                    if category in expected_categories:
                        expected_energy += find_value(category)
                    elif category_energy:
                        raise ValueError(f"expected zero energy from CHARMM for {category}, not {category_energy}")
            if abs(expected_energy - energy) > CHARMM_ENERGY_CHECK_PRECISION:
                print(result.stdout.decode())
                raise ValueError(
                    f"expected total energy {energy:.5f} from CHARMM to match sum {expected_energy:.5f} of terms"
                )

            # Look for the forces.
            advance_until(lambda line: "COORDINATE" in line)
            advance_until(lambda line: "TITLE>" not in line)

            # Extract the forces.
            forces = []
            for atom_index in range(len(coordinates_list[coordinates_index])):
                force_data = next(lines).split()
                if int(force_data[0]) != atom_index + 1:
                    raise ValueError(f"unexpected force data at atom {atom_index}")
                forces.append(list(map(float, force_data[4:7])))

            return dict(
                energy=energy * openmm.unit.kilocalorie_per_mole,
                # The CHARMM forces are inverted when read, so flip them.
                forces=openmm.unit.Quantity(
                    -numpy.array(forces), openmm.unit.kilocalorie_per_mole / openmm.unit.angstrom
                ),
                mask=virtual_site_mask,
            )

        results = []
        for coordinates_index in range(len(coordinates_list)):
            result = {force_group: parse_energies_forces(force_group) for force_group in (None, *ForceGroup)}

            # Rescale the electrostatic energy and forces by the ratio of the vacuum
            # permittivities used in CHARMM and OpenMM.
            electrostatic_charmm = parse_energies_forces(None)
            old_energy = electrostatic_charmm["energy"]
            old_forces = electrostatic_charmm["forces"]
            new_energy = old_energy * CHARMM_ELECTROSTATIC_SCALE
            new_forces = old_forces * CHARMM_ELECTROSTATIC_SCALE

            for force_group in (None, ForceGroup.NONBONDED):
                result[force_group]["energy"] += new_energy - old_energy
                result[force_group]["forces"] += new_forces - old_forces

            results.append(result)

        return results

    def get_openmm_charmm_energies_forces(self, test_directory, test_spec, coordinates_list, dump_path):
        """
        Computes energies and forces with OpenMM using OpenMM CharmmParameterSet.

        Parameters
        ----------
        test_directory : str
            Path to the directory containing a test specification YAML file.
        test_spec : dict
            A test specification read from a YAML file.
        coordinates_list : list
            A list of arrays of positions.  Energies and forces will be evaluated
            for each set of coordinates.
        dump_path : str or None
            Path to dump input files for debugging, or None to skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary of dictionaries for each
            ForceGroup (or None for all force field terms) containing the system
            energy as "energy" and forces as "forces".
        """

        _, system = self._generate_psf_system(test_directory, test_spec)
        return self._evaluate_openmm(system, coordinates_list, dump_path, test_spec.get("openmm_pme", None))

    def get_parmed_charmm_energies_forces(self, test_directory, test_spec, coordinates_list, dump_path):
        """
        Computes energies and forces with OpenMM using ParmEd CharmmParameterSet.

        Parameters
        ----------
        test_directory : str
            Path to the directory containing a test specification YAML file.
        test_spec : dict
            A test specification read from a YAML file.
        coordinates_list : list
            A list of arrays of positions.  Energies and forces will be evaluated
            for each set of coordinates.
        dump_path : str or None
            Path to dump input files for debugging, or None to skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary of dictionaries for each
            ForceGroup (or None for all force field terms) containing the system
            energy as "energy" and forces as "forces".
        """

        # ParmEd is not relied on anywhere else, so this test can be skipped
        # safely if ParmEd is not installed.
        import parmed

        parameter_set = parmed.charmm.CharmmParameterSet(*self._get_toppar_files(test_spec))
        psf_file = parmed.charmm.CharmmPsfFile(os.path.join(test_directory, test_spec["psf_file"]))

        box = self._get_openmm_box(test_spec)
        if box is not None:
            psf_file.box = (box[0, 0], box[1, 1], box[2, 2], 90, 90, 90)

        system = psf_file.createSystem(parameter_set, **test_spec["create_system_options"])
        return self._evaluate_openmm(system, coordinates_list, dump_path, test_spec.get("openmm_pme", None))

    def get_openmm_ffxml_energies_forces(self, test_directory, test_spec, coordinates_list, dump_path):
        """
        Computes energies and forces with OpenMM using OpenMM ForceField.

        Parameters
        ----------
        test_directory : str
            Path to the directory containing a test specification YAML file.
        test_spec : dict
            A test specification read from a YAML file.
        coordinates_list : list
            A list of arrays of positions.  Energies and forces will be evaluated
            for each set of coordinates.
        dump_path : str or None
            Path to dump input files for debugging, or None to skip dumping.

        Returns
        -------
        list
            For each set of coordinates, a dictionary of dictionaries for each
            ForceGroup (or None for all force field terms) containing the system
            energy as "energy" and forces as "forces".
        """

        force_field = openmm.app.ForceField(
            *(os.path.join(self.ffxml_directory, file) for file in test_spec["ffxml_files"])
        )
        pdb_file = openmm.app.PDBFile(os.path.join(test_directory, test_spec["pdb_file"]))
        topology = pdb_file.getTopology()

        box = self._get_openmm_box(test_spec)
        if box is not None:
            topology.setPeriodicBoxVectors(box)

        extra_options = {}
        if "ffxml_dispersion_correction" in test_spec:
            extra_options["useDispersionCorrection"] = test_spec["ffxml_dispersion_correction"]
        system = force_field.createSystem(topology, **test_spec["create_system_options"], **extra_options)

        if self.openmm_ffxml_fix_impropers:
            psf_file, psf_system = self._generate_psf_system(test_directory, test_spec)
            self._fix_impropers(psf_file, psf_system, system)

        return self._evaluate_openmm(system, coordinates_list, dump_path, test_spec.get("openmm_pme", None))

    def _generate_psf_system(self, test_directory, test_spec):
        """
        Uses OpenMM to generate a System from a CHARMM PSF file.

        Parameters
        ----------
        test_directory : str
            Path to the directory containing a test specification YAML file.
        test_spec : dict
            A test specification read from a YAML file.

        Returns
        -------
        openmm.app.CharmmPsfFile, openmm.System
            The PSF file loaded, and the OpenMM system created.
        """

        parameter_set = openmm.app.CharmmParameterSet(*self._get_toppar_files(test_spec))
        psf_file = openmm.app.CharmmPsfFile(
            os.path.join(test_directory, test_spec["psf_file"]), periodicBoxVectors=self._get_openmm_box(test_spec)
        )
        return psf_file, psf_file.createSystem(parameter_set, **test_spec["create_system_options"])

    def _fix_impropers(self, psf_file, psf_system, ffxml_system):
        """
        Tries to correct the order of the impropers in an OpenMM system to be
        consistent with a PSF file.  Makes sure that the OpenMM FFXML impropers
        are as consistent as possible with the PSF impropers.

        Parameters
        ----------
        psf_file : openmm.app.CharmmPsfFile
            A PSF file containing reference data about the system's impropers.
        psf_system : openmm.System
            A system created from the PSF file.
        ffxml_system : openmm.System
            A system created from an OpenMM FFXML force field.
        """

        # Determine which atoms are bonded to which other atoms, to determine
        # which torsions are improper torsions.
        atoms_bonded_to = [set() for index in range(len(psf_file.atom_list))]
        for bond in psf_file.bond_list:
            atoms_bonded_to[bond.atom1.idx].add(bond.atom2.idx)
            atoms_bonded_to[bond.atom2.idx].add(bond.atom1.idx)

        # Helper function to find if atom indices correspond to an improper.
        # The central atom should always be listed first or last.
        def check_improper(index_1, index_2, index_3, index_4):
            bonded_to_1 = atoms_bonded_to[index_1]
            bonded_to_4 = atoms_bonded_to[index_4]
            is_improper_1 = index_2 in bonded_to_1 and index_3 in bonded_to_1 and index_4 in bonded_to_1
            is_improper_4 = index_1 in bonded_to_4 and index_2 in bonded_to_4 and index_3 in bonded_to_4

            # Return whether or not this is an improper, and if so, if the
            # central atom is last instead of first.  If both atoms look like a
            # central atom, something is wrong; report that this is not an
            # improper.  If this is not an improper, the second return value
            # should not be checked.
            return (is_improper_1 ^ is_improper_4, is_improper_4)

        # Helper function to turn a set of indices into a key for an improper.
        def improper_to_key(index_1, index_2, index_3, index_4, is_flipped):
            if is_flipped:
                return (index_4, *sorted((index_1, index_2, index_3)))
            else:
                return (index_1, *sorted((index_2, index_3, index_4)))

        # Classify the impropers in the PSF.
        psf_impropers = {}
        for improper in psf_file.improper_list:
            index_1 = improper.atom1.idx
            index_2 = improper.atom2.idx
            index_3 = improper.atom3.idx
            index_4 = improper.atom4.idx

            # Make sure that this improper is topologically valid.
            is_improper, is_flipped = check_improper(index_1, index_2, index_3, index_4)
            if not is_improper:
                raise ValueError(f"improper {index_1}-{index_2}-{index_3}-{index_4} in PSF is invalid")

            # Make a key for the improper lookup that ignores the order of the
            # non-central atoms.  Store in a dictionary the correct order.
            # There may be more than one improper for this key.
            psf_impropers.setdefault(improper_to_key(index_1, index_2, index_3, index_4, is_flipped), []).append(
                (index_1, index_2, index_3, index_4)
            )

        # Classifies the impropers in a system.
        def classify_system_impropers(system):
            impropers = {}
            for force_index, force in enumerate(system.getForces()):
                if isinstance(force, (openmm.PeriodicTorsionForce, openmm.CustomTorsionForce)):
                    for term_index in range(force.getNumTorsions()):
                        index_1, index_2, index_3, index_4, *torsion_parameters = force.getTorsionParameters(
                            term_index
                        )
                        is_improper, is_flipped = check_improper(index_1, index_2, index_3, index_4)
                        if not is_improper:
                            continue
                        impropers.setdefault(
                            improper_to_key(index_1, index_2, index_3, index_4, is_flipped), []
                        ).append(((index_1, index_2, index_3, index_4), (force_index, term_index)))
            return impropers

        psf_system_impropers = classify_system_impropers(psf_system)
        ffxml_system_impropers = classify_system_impropers(ffxml_system)

        # Make sure that the PSF system impropers match the PSF impropers
        # perfectly.
        if set(psf_impropers.keys()) != set(psf_system_impropers.keys()) or any(
            sorted(psf_impropers[key])
            != sorted(atom_indices for atom_indices, system_indices in psf_system_impropers[key])
            for key in psf_impropers.keys()
        ):
            raise ValueError("PSF and PSF system impropers do not match")

        # Make sure that the FFXML system impropers match as closely as
        # possible.  We ensure that the correct sets of atoms and numbers of
        # permutations are present.
        if set(psf_system_impropers.keys()) != set(ffxml_system_impropers.keys()) or any(
            len(psf_system_impropers[key]) != len(ffxml_system_impropers[key]) for key in psf_system_impropers.keys()
        ):
            raise ValueError("PSF and FFXML system impropers do not match")

        for key, psf_list in psf_system_impropers.items():
            improper_count = len(psf_list)
            ffxml_list = ffxml_system_impropers.get(key, [])

            if improper_count == 1:
                # There is a single improper for the central atom and set of
                # non-central atoms.
                ((psf_atom_indices, (psf_force_index, psf_term_index)),) = psf_list
                ((ffxml_atom_indices, (ffxml_force_index, ffxml_term_index)),) = ffxml_list

                # See if there is nothing wrong with this improper.
                permuted_indices = [
                    psf_atom_index
                    for psf_atom_index, ffxml_atom_index in zip(psf_atom_indices, ffxml_atom_indices)
                    if psf_atom_index != ffxml_atom_index
                ]
                if len(permuted_indices) < 2:
                    continue

                # Make sure that if atoms are permuted, the permuted atoms
                # appear chemically equivalent (same classes, and bonded to
                # atoms of the same classes).
                reference_index, *permuted_indices = permuted_indices
                reference_types = sorted(
                    psf_file.atom_list[bonded_index].attype for bonded_index in atoms_bonded_to[reference_index]
                )
                for permuted_index in permuted_indices:
                    if psf_file.atom_list[reference_index].attype != psf_file.atom_list[permuted_index].attype:
                        raise ValueError("Permuted atoms in FFXML system have inequivalent classes")
                    permuted_types = sorted(
                        psf_file.atom_list[bonded_index].attype for bonded_index in atoms_bonded_to[permuted_index]
                    )
                    if reference_types != permuted_types:
                        raise ValueError(
                            "Permuted atoms in FFXML system are bonded to atoms with inequivalent classes"
                        )

                # If this point is reached, it is safe to try to replace the
                # permuted FFXML improper with the indices from the PSF system.
                psf_torsion_parameters = psf_system.getForce(psf_force_index).getTorsionParameters(psf_term_index)[4:]
                ffxml_torsion_parameters = ffxml_system.getForce(ffxml_force_index).getTorsionParameters(
                    ffxml_term_index
                )[4:]
                if psf_torsion_parameters != ffxml_torsion_parameters:
                    raise ValueError("PSF and FFXML system improper parameters do not match")
                ffxml_system.getForce(ffxml_force_index).setTorsionParameters(
                    ffxml_term_index, *psf_atom_indices, *ffxml_torsion_parameters
                )

            else:
                # If there is more than one improper for a central atom and set
                # of non-central atoms, give up unless there is a perfect match.
                if not all(
                    psf_atom_indices == ffxml_atom_indices
                    for (psf_atom_indices, _), (ffxml_atom_indices, _) in zip(psf_list, ffxml_list)
                ):
                    raise ValueError("Don't know how to make PSF and FFXML system impropers correspond")

    def _get_toppar_files(self, test_spec):
        """
        Extracts a list of all topology, parameter, and stream files from a test
        specification.

        Parameters
        ----------
        test_spec : dict
            A test specification read from a YAML file.

        Returns
        -------
        list
            A list of paths to files in the CHARMM parameter set directory.
        """

        return [os.path.join(self.toppar_directory, file) for file in test_spec["charmm_files"]]

    def _get_openmm_box(self, test_spec):
        """
        Converts a box specification to a set of OpenMM box vectors.

        Parameters
        ----------
        test_spec : dict
            A test specification read from a YAML file.

        Returns
        -------
        openmm.unit.Quantity or None
            Orthogonal box vectors as a diagonal matrix, or None if the test
            does not take place in a periodic system.
        """

        if "box" not in test_spec:
            return None

        box = test_spec["box"]
        return (
            numpy.diag(
                (
                    box["a"].value_in_unit(openmm.unit.angstrom),
                    box["b"].value_in_unit(openmm.unit.angstrom),
                    box["c"].value_in_unit(openmm.unit.angstrom),
                )
            )
            * openmm.unit.angstrom
        )

    def _evaluate_openmm(self, system, coordinates_list, dump_path, pme_parameters):
        """
        Evaluates energies and forces for an OpenMM system.  This will modify
        the system.

        Parameters
        ----------
        system : openmm.System
            The system for which to create a context and perform the evaluation.
        coordinates_list : list
            A list of arrays of positions.  Energies and forces will be
            evaluated for each set of coordinates.
        dump_path : str or None
            Path to dump input files for debugging, or None to skip dumping.
        pme_parameters: dict or None
            Parameters fftx, ffty, fftz, and kappa used to match the PME
            behavior between CHARMM and OpenMM.  If None, OpenMM will be allowed
            to choose its own PME parameters.

        Returns
        -------
        list
            For each set of coordinates, a dictionary of dictionaries for each
            ForceGroup (or None for all force field terms) containing the system
            energy as "energy" and forces as "forces".
        """

        # Get a mask of virtual sites.
        virtual_site_mask = numpy.array(
            [system.isVirtualSite(particle_index) for particle_index in range(system.getNumParticles())], dtype=bool
        )

        # Set force groups of all forces in the System.
        for force_index, force in enumerate(system.getForces()):
            for force_group in ForceGroup:
                if isinstance(force, FORCE_GROUP_DATA[force_group][2]):
                    force.setForceGroup(force_group.value)
                    break
            else:
                raise RuntimeError(f"unknown force {type(force).__name__}")

        if pme_parameters is not None:
            for force in system.getForces():
                if isinstance(force, openmm.NonbondedForce):
                    force.setPMEParameters(
                        1 / pme_parameters["kappa"],
                        pme_parameters["fftx"],
                        pme_parameters["ffty"],
                        pme_parameters["fftz"],
                    )

        if dump_path is not None:
            with open(dump_path, "w") as dump_file:
                dump_file.write(openmm.XmlSerializer.serialize(system))

        # Dummy integrator should not be utilized.
        integrator = openmm.VerletIntegrator(0.001)

        if self.openmm_platform == "automatic":
            context = openmm.Context(system, integrator)
        else:
            context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName(self.openmm_platform))

        results = []
        for coordinates in coordinates_list:
            context.setPositions(coordinates)
            context.computeVirtualSites()

            result = {}

            # Evaluate energies and forces for the entire system.
            state = context.getState(getEnergy=True, getForces=True)
            result[None] = dict(
                energy=state.getPotentialEnergy(), forces=state.getForces(asNumpy=True), mask=virtual_site_mask
            )

            # Evaluate energies and forces for each force group.
            for force_group in ForceGroup:
                state = context.getState(getEnergy=True, getForces=True, groups={force_group.value})
                result[force_group] = dict(
                    energy=state.getPotentialEnergy(), forces=state.getForces(asNumpy=True), mask=virtual_site_mask
                )

            results.append(result)

        return results


def random_in_unit_ball(generator):
    """
    Generates a random point in the closed unit ball.

    Parameters
    ----------
    generator : numpy.random.Generator
        Generator used to produce random values.

    Returns
    -------
    numpy.ndarray
        3-vector of floating-point numbers uniformly randomly drawn from the
        closed unit ball.
    """

    # Generate by rejection sampling starting from an axis-aligned cube of side
    # length 2 centered at the origin.
    while True:
        vector = generator.uniform(-1, 1, 3)
        if vector @ vector <= 1:
            return vector


def color_message(message, ansi="1;38;5;216"):
    """
    Generates ANSI escape sequences to highlight a message with color or other
    effects.

    Parameters
    ----------
    message : str
        The message to color.
    ansi : str
        The ANSI formatting codes to apply.

    Returns
    -------
    str
        A colored message, if the standard output stream is a TTY; otherwise,
        the message unmodified.
    """

    return f"\x1b[{ansi}m{message}\x1b[0m" if DO_COLOR else message


def yaml_to_quantity(yaml_quantity):
    """
    Converts a quantity as [value, unit_name] read from a YAML file into an
    OpenMM Quantity.

    Parameters
    ----------
    yaml_quantity : (float, str)
        A specification of a quantity read from a YAML file.

    Returns
    -------
    openmm.unit.Quantity
        A corresponding quantity.
    """

    value, unit_name = yaml_quantity
    unit = None
    if hasattr(openmm.unit.unit_definitions, unit_name):
        unit_definition = getattr(openmm.unit.unit_definitions, unit_name)
        if openmm.unit.is_unit(unit_definition):
            unit = unit_definition
    if unit is None:
        raise ValueError(f"unknown unit {unit_name!r}")
    return value * unit


if __name__ == "__main__":
    main()
