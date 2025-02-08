import argparse
import enum
import itertools
import numpy
import openmm.app
import os
import subprocess
import sys
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

DEFAULT_ABSOLUTE_ENERGY_TOLERANCE = 5e-5  # kcal/mol
DEFAULT_RELATIVE_ENERGY_TOLERANCE = 5e-5
DEFAULT_ABSOLUTE_FORCE_TOLERANCE = 5e-3  # kcal/mol/Å
DEFAULT_RELATIVE_FORCE_TOLERANCE = 5e-3

DEFAULT_PERTURB_DISTANCE = 0.1  # Å
DEFAULT_TOPPAR_DIRECTORY = "toppar"
DEFAULT_FFXML_DIRECTORY = "ffxml"

ENERGY_PRINT_UNIT = openmm.unit.kilocalorie_per_mole
FORCE_PRINT_UNIT = ENERGY_PRINT_UNIT / openmm.unit.angstrom

# A string that shouldn't appear in the CHARMM output.
ENERGY_FORCE_TAG = "OPENMMFORCEFIELDS-TEST-CHARMM-ENERGY-FORCES"

# CHARMM doesn't print out energies to more than 5 places, so the last place
# might be off by a few digits.
CHARMM_ENERGY_CHECK_PRECISION = 5e-5

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


# For each force group, the CHARMM skip keywords, the labels in the CHARMM
# output, and the OpenMM force types.
FORCE_GROUP_DATA = {
    ForceGroup.BONDS_UREY_BRADLEY: (("bond", "urey"), ("BONDs", "UREY-b"), (openmm.HarmonicBondForce,)),
    ForceGroup.ANGLES: (("angl",), ("ANGLes",), (openmm.HarmonicAngleForce,)),
    ForceGroup.PROPERS: (("dihe",), ("DIHEdrals",), (openmm.PeriodicTorsionForce,)),
    ForceGroup.IMPROPERS: (("impr",), ("IMPRopers",), (openmm.CustomTorsionForce,)),
    ForceGroup.CMAP: (("cmap",), ("CMAPs",), (openmm.CMAPTorsionForce,)),
    ForceGroup.NONBONDED: (
        ("vdw", "elec"),
        ("VDWaals", "ELEC"),
        (openmm.NonbondedForce, openmm.CustomNonbondedForce, openmm.CustomBondForce),
    ),
    # NOTE: Impropers could be under a PeriodicTorsionForce but this almost
    # never occurs in the CHARMM force field.  The decomposition between propers
    # and impropers would fail in that case even if the total energy was being
    # evaluated correctly, but this shouldn't happen for most systems.
}


def main():
    """
    Main entry point for test script.  Parses command-line arguments and runs.
    """

    parser = argparse.ArgumentParser(
        description="Test script for CHARMM force field conversion to OpenMM FFXML format"
    )

    parser.add_argument(
        "test_directory",
        nargs="*",
        help="Directory in which to find a test specification (defined by a test.yaml file)",
    )

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
        "--openmm-ffxml", action="store_true", help="Compute energies and forces with OpenMM using OpenMM ForceField"
    )

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
        "--perturb-seed", type=int, metavar="seed", help="Random seed for repeating tests with perturbed positions"
    )

    parser.add_argument(
        "--toppar-directory", default="toppar", metavar="path", help="Path to CHARMM parameter files (default toppar)"
    )
    parser.add_argument(
        "--ffxml-directory", default="ffxml", metavar="path", help="Path to OpenMM force field files (default ffxml)"
    )

    parser.add_argument("--dump", action="store_true", help="Dump input files for debugging")

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

    # Run all tests, storing (test_directory, success_flag, message) for each.
    test_results = []
    for test_directory in args.test_directory:
        try:
            term_failure_count = run_test(
                args.toppar_directory,
                args.ffxml_directory,
                test_directory,
                args.charmm,
                args.openmm_charmm,
                args.parmed_charmm,
                args.openmm_ffxml,
                perturb_spec,
                args.absolute_energy_tolerance * openmm.unit.kilocalorie_per_mole,
                args.relative_energy_tolerance,
                args.absolute_force_tolerance * openmm.unit.kilocalorie_per_mole / openmm.unit.angstrom,
                args.relative_force_tolerance,
                skip_set,
                args.dump,
            )
            test_results.append(
                (test_directory, term_failure_count == 0, f"{term_failure_count} terms exceeded tolerance")
            )
        except Exception as error:
            test_results.append((test_directory, False, f"{type(error).__name__}: {error}"))

    if not test_results:
        print("No tests run")
        sys.exit()

    # Summarize test results.
    test_count = len(test_results)
    test_names, success_flags, _ = zip(*test_results)
    success_count = sum(success_flags)
    failure_count = test_count - success_count

    print(f"Ran {test_count} tests ({success_count} successes, {failure_count} failures)")
    test_name_width = max(map(len, test_names))
    for test_name, success_flag, test_message in test_results:
        formatted_line = f"{"Succeeded" if success_flag else "Failed":9}  {test_message}"
        print(
            f"    {test_name:{test_name_width}}  {formatted_line if success_flag else color_message(formatted_line)}"
        )

    sys.exit(failure_count)


def run_test(
    toppar_directory,
    ffxml_directory,
    test_directory,
    do_charmm,
    do_openmm_charmm,
    do_parmed_charmm,
    do_openmm_ffxml,
    perturb_spec,
    absolute_energy_tolerance,
    relative_energy_tolerance,
    absolute_force_tolerance,
    relative_force_tolerance,
    skip_set,
    dump,
):
    """
    Runs the specified test with the specified options.

    Parameters
    ----------
    toppar_directory : str
        Directory containing CHARMM parameter files.
    ffxml_directory : str
        Directory containing OpenMM force field files.
    test_directory : str
        Directory in which to find a test specification (defined by a test.yaml
        file).
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
    perturb_spec : (numpy.random.Generator, int, float) or None
        If not None, repeat tests the given number of times with positions
        perturbed by up to the given distance in angstroms using random numbers
        from the given generator.
    absolute_energy_tolerance : openmm.unit.Quantity
        Absolute tolerance for matching energies.
    relative_energy_tolerance : float
        Relative tolerance for matching energies.
    absolute_force_tolerance : openmm.unit.Quantity
        Absolute tolerance for matching forces.
    relative_force_tolerance : float
        Relative tolerance for matching forces.
    skip_set : set(ForceGroup)
        Set of terms to ignore during testing.
    dump : bool
        Whether or not to dump input files for debugging.

    Returns
    -------
    int
        The number of failed checks.
    """

    with open(os.path.join(test_directory, "test.yaml")) as test_spec_file:
        test_spec = yaml.safe_load(test_spec_file)

    # Generate reference coordinates.
    coordinates_list = [openmm.app.PDBFile(os.path.join(test_directory, test_spec["pdb_file"])).getPositions()]
    if perturb_spec is not None:
        generator, perturb_replicates, perturb_distance = perturb_spec
        for replicate in range(perturb_replicates):
            perturbed_coordinates = coordinates_list[0].value_in_unit(openmm.unit.angstrom).copy()
            for index in range(len(perturbed_coordinates)):
                perturbed_coordinates[index] += openmm.Vec3(*(perturb_distance * random_in_unit_ball(generator)))
            coordinates_list.append(openmm.unit.Quantity(perturbed_coordinates, openmm.unit.angstrom))

    dump_path = test_directory.replace("/", "_")

    # TODO: need to get and pass periodic box information from the spec file to
    # CHARMM and OpenMM.

    # Run tests and get results.
    results_sets = {}
    if do_charmm:
        results_sets["charmm"] = get_charmm_energies_forces(
            toppar_directory, test_directory, test_spec, coordinates_list, f"{dump_path}_charmm.inp" if dump else None
        )
    if do_openmm_charmm:
        results_sets["openmm_charmm"] = get_openmm_charmm_energies_forces(
            toppar_directory,
            test_directory,
            test_spec,
            coordinates_list,
            f"{dump_path}_openmm_charmm.xml" if dump else None,
        )
    if do_parmed_charmm:
        results_sets["parmed_charmm"] = get_parmed_charmm_energies_forces(
            toppar_directory,
            test_directory,
            test_spec,
            coordinates_list,
            f"{dump_path}_parmed_charmm.xml" if dump else None,
        )
    if do_openmm_ffxml:
        results_sets["openmm_ffxml"] = get_openmm_ffxml_energies_forces(
            ffxml_directory,
            test_directory,
            test_spec,
            coordinates_list,
            f"{dump_path}_openmm_ffxml.xml" if dump else None,
        )

    # Zero anything that's not being tested.
    for results in results_sets.values():
        for result in results:
            for skip_group in skip_set:
                result[None]["energy"] -= result[skip_group]["energy"]
                result[skip_group]["energy"] *= 0
                result[None]["forces"] -= result[skip_group]["forces"]
                result[skip_group]["forces"] *= 0

    # Helper function for comparing energies.
    def compare_energies(energy_1, energy_2):
        absolute_difference = abs(energy_2 - energy_1)
        if energy_1 or energy_2:
            relative_difference = absolute_difference / max(abs(energy_1), abs(energy_2))
        else:
            relative_difference = 0.0
        return absolute_difference, relative_difference

    # Helper function for comparing forces.
    def compare_forces(force_1, force_2):
        magnitude_1 = numpy.linalg.norm(force_1.value_in_unit(FORCE_PRINT_UNIT)) * FORCE_PRINT_UNIT
        magnitude_2 = numpy.linalg.norm(force_2.value_in_unit(FORCE_PRINT_UNIT)) * FORCE_PRINT_UNIT
        absolute_difference = numpy.linalg.norm((force_2 - force_1).value_in_unit(FORCE_PRINT_UNIT)) * FORCE_PRINT_UNIT
        if magnitude_1 or magnitude_2:
            relative_difference = absolute_difference / max(magnitude_1, magnitude_2)
        else:
            relative_difference = 0.0
        return (
            absolute_difference > absolute_force_tolerance and relative_difference > relative_force_tolerance,
            absolute_difference,
            relative_difference,
        )

    # Helper function for printing values.
    def format_energy_difference_data(label, difference_data):
        energies_1, energies_2, absolute_differences, relative_differences = zip(*difference_data)
        failures = [
            absolute_difference > absolute_energy_tolerance and relative_difference > relative_energy_tolerance
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
                print(f"{" " * 66}{color_message(formatted_line) if failure else formatted_line}")

        return any_failures

    def format_force_difference_data(label, difference_data):
        failures, absolute_differences, relative_differences = zip(*difference_data)

        any_failures = any(failures)
        formatted_line = (
            f"{max(absolute_differences).value_in_unit(FORCE_PRINT_UNIT):20.12e}  {max(relative_differences):12.8f}"
        )
        if any_failures:
            formatted_line = color_message(formatted_line)
        print(f"        {label:20}  {formatted_line}")

        return any_failures

    if len(results_sets) < 2:
        raise ValueError("must specify at least 2 methods to compare")

    print(f"Ran test {test_spec["name"]!r} at {test_directory}")
    failure_count = 0
    for (name_1, results_1), (name_2, results_2) in itertools.combinations(results_sets.items(), 2):
        print(f"    Comparing {name_1} vs. {name_2}")

        # Print results for energies.
        print(f"        {"Energy error":20}  {"|ΔE| (kcal/mol)":20}  {"Relative":12}")
        for force_group in (None, *ForceGroup):
            difference_data = []
            for result_1, result_2 in zip(results_1, results_2):
                energy_1 = result_1[force_group]["energy"]
                energy_2 = result_2[force_group]["energy"]
                difference_data.append((energy_1, energy_2, *compare_energies(energy_1, energy_2)))
            failure_count += format_energy_difference_data(
                "TOTAL" if force_group is None else force_group.name, difference_data
            )
        print()

        # Print results for forces.
        print(f"        {"Force error":20}  {"|Δf| (kcal/mol/Å)":20}  {"Relative":12}")
        for force_group in (None, *ForceGroup):
            difference_data = []
            for result_1, result_2 in zip(results_1, results_2):
                for force_1, force_2 in zip(result_1[force_group]["forces"], result_2[force_group]["forces"]):
                    difference_data.append(compare_forces(force_1, force_2))
            failure_count += format_force_difference_data(
                "TOTAL" if force_group is None else force_group.name, difference_data
            )
        print()

    return failure_count


def get_charmm_energies_forces(toppar_directory, test_directory, test_spec, coordinates_list, dump_path):
    """
    Computes energies and forces with CHARMM.

    Parameters
    ----------
    toppar_directory : str
        Directory containing CHARMM parameter files.
    test_directory : str
        Path to the directory containing a test.yaml file.
    test_spec : dict
        A test specification read from a test.yaml file.
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

    # Prepare the CHARMM input line by line.
    charmm_input_lines = ["* Test", "*"]

    # Load all input files.
    for rtf_file in test_spec["rtf_files"]:
        charmm_input_lines.append(f"read rtf card name {os.path.join(toppar_directory, rtf_file)}")
    for prm_file in test_spec["prm_files"]:
        charmm_input_lines.append(f"read para card name {os.path.join(toppar_directory, prm_file)} flex")
    for str_file in test_spec["str_files"]:
        charmm_input_lines.append(f"stream {os.path.join(toppar_directory, str_file)}")
    charmm_input_lines.append(f"read psf card name {os.path.join(test_directory, test_spec["psf_file"])}")

    # Execute custom commands.
    charmm_input_lines.extend(test_spec["charmm_commands"])

    for coordinates_index, coordinates in enumerate(coordinates_list):
        # Set all coordinates.
        for atom_index, coordinate in enumerate(coordinates):
            x, y, z = coordinate.value_in_unit(openmm.unit.angstrom)
            charmm_input_lines.append(f"coor set xdir {x} ydir {y} zdir {z} sele bynu {atom_index + 1} end")

        # Evaluate the energy and forces for the entire system.  Append a tag
        # that we can search for in the output.
        charmm_input_lines.append(f"!{ENERGY_FORCE_TAG}")
        charmm_input_lines.append("skip none")
        charmm_input_lines.append("ener")
        charmm_input_lines.append("coor force comp")
        charmm_input_lines.append("print coor comp")

        # Evaluate energies and forces for each force group.  Turn off all terms
        # (skip all), but leave (excl) some in.
        for force_group in ForceGroup:
            charmm_input_lines.append(f"! {ENERGY_FORCE_TAG}")
            charmm_input_lines.append(f"skip all excl {" ".join(FORCE_GROUP_DATA[force_group][0])}")
            charmm_input_lines.append("ener")
            charmm_input_lines.append("coor force comp")
            charmm_input_lines.append("print coor comp")

        # Evaluate just electrostatic energies and forces.  This is used to get
        # a good comparison between OpenMM and CHARMM because of their different
        # values for the vacuum permittivity.
        charmm_input_lines.append(f"! {ENERGY_FORCE_TAG}")
        charmm_input_lines.append("skip all excl elec")
        charmm_input_lines.append("ener")
        charmm_input_lines.append("coor force comp")
        charmm_input_lines.append("print coor comp")

    if dump_path is not None:
        with open(dump_path, "w") as dump_file:
            for line in charmm_input_lines:
                print(line, file=dump_file)

    # Run CHARMM.
    charmm_input_lines.append("stop")
    result = subprocess.run(["charmm"], input="\n".join(charmm_input_lines).encode(), capture_output=True)

    # If something went wrong, print the output for debugging.
    if result.returncode:
        print(result.stdout.decode())
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
        ener_values = [values.strip().split() for values in ener_values]

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

            # The first line of labels starts with "ENER ENR:" and the first
            # line of values starts with "ENER>", so this has to be adjusted.
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
            forces=openmm.unit.Quantity(-numpy.array(forces), openmm.unit.kilocalorie_per_mole / openmm.unit.angstrom),
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


def get_openmm_charmm_energies_forces(toppar_directory, test_directory, test_spec, coordinates_list, dump_path):
    """
    Computes energies and forces with OpenMM using OpenMM CharmmParameterSet.

    Parameters
    ----------
    toppar_directory : str
        Directory containing CHARMM parameter files.
    test_directory : str
        Path to the directory containing a test.yaml file.
    test_spec : dict
        A test specification read from a test.yaml file.
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

    parameter_set = openmm.app.CharmmParameterSet(*get_toppar_files(toppar_directory, test_spec))
    psf_file = openmm.app.CharmmPsfFile(os.path.join(test_directory, test_spec["psf_file"]))
    system = psf_file.createSystem(parameter_set, **test_spec["create_system_options"])
    return evaluate_openmm(system, coordinates_list, dump_path)


def get_parmed_charmm_energies_forces(toppar_directory, test_directory, test_spec, coordinates_list, dump_path):
    """
    Computes energies and forces with OpenMM using ParmEd CharmmParameterSet.

    Parameters
    ----------
    toppar_directory : str
        Directory containing CHARMM parameter files.
    test_directory : str
        Path to the directory containing a test.yaml file.
    test_spec : dict
        A test specification read from a test.yaml file.
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

    # ParmEd is not relied on anywhere else, so this test can be skipped safely
    # if ParmEd is not installed.
    import parmed

    parameter_set = parmed.charmm.CharmmParameterSet(*get_toppar_files(toppar_directory, test_spec))
    psf_file = parmed.charmm.CharmmPsfFile(os.path.join(test_directory, test_spec["psf_file"]))
    system = psf_file.createSystem(parameter_set, **test_spec["create_system_options"])
    return evaluate_openmm(system, coordinates_list, dump_path)


def get_openmm_ffxml_energies_forces(ffxml_directory, test_directory, test_spec, coordinates_list, dump_path):
    """
    Computes energies and forces with OpenMM using OpenMM ForceField.

    Parameters
    ----------
    ffxml_directory : str
        Directory containing CHARMM parameter files.
    test_directory : str
        Path to the directory containing a test.yaml file.
    test_spec : dict
        A test specification read from a test.yaml file.
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

    force_field = openmm.app.ForceField(*(os.path.join(ffxml_directory, file) for file in test_spec["ffxml_files"]))
    pdb_file = openmm.app.PDBFile(os.path.join(test_directory, test_spec["pdb_file"]))
    system = force_field.createSystem(pdb_file.getTopology(), **test_spec["create_system_options"])
    return evaluate_openmm(system, coordinates_list, dump_path)


def get_toppar_files(toppar_directory, test_spec):
    """
    Extracts a list of all topology, parameter, and stream files from a test
    specification.  This is done for convenience since OpenMM and ParmEd do not
    require separate handling for the different kinds of files.

    Parameters
    ----------
    toppar_directory : str
        Directory containing CHARMM parameter files.
    test_spec : dict
        A test specification read from a test.yaml file.

    Returns
    -------
    list
        A list of paths to files in the CHARMM parameter set directory.
    """

    return [
        os.path.join(toppar_directory, file)
        for file in (*test_spec["rtf_files"], *test_spec["prm_files"], *test_spec["str_files"])
    ]


def evaluate_openmm(system, coordinates_list, dump_path):
    """
    Evaluates energies and forces for an OpenMM system.  This will modify the
    system.

    Parameters
    ----------
    system : openmm.System
        The system for which to create a context and perform the evaluation.
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

    # Set force groups of all forces in the System.
    for force_index, force in enumerate(system.getForces()):
        for force_group in ForceGroup:
            if isinstance(force, FORCE_GROUP_DATA[force_group][2]):
                force.setForceGroup(force_group.value)
                break
        else:
            raise RuntimeError(f"unknown force {type(force).__name__}")

    if dump_path is not None:
        with open(dump_path, "w") as dump_file:
            dump_file.write(openmm.XmlSerializer.serialize(system))

    context = openmm.Context(system, openmm.VerletIntegrator(0.001))

    results = []
    for coordinates in coordinates_list:
        context.setPositions(coordinates)

        result = {}

        # Evaluate energies and forces for the entire system.
        state = context.getState(getEnergy=True, getForces=True)
        result[None] = dict(energy=state.getPotentialEnergy(), forces=state.getForces(asNumpy=True))

        # Evaluate energies and forces for each force group.
        for force_group in ForceGroup:
            state = context.getState(getEnergy=True, getForces=True, groups={force_group.value})
            result[force_group] = dict(energy=state.getPotentialEnergy(), forces=state.getForces(asNumpy=True))

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


if __name__ == "__main__":
    main()
