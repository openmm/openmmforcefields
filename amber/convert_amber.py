# AMBER --> OpenMM force-field conversion script
# Author: Rafal P. Wiewiora, ChoderaLab
import argparse
import csv
import glob
import hashlib
import os
import re
import sys
import tempfile
import warnings
import xml.etree.ElementTree as etree
from copy import deepcopy
from distutils.spawn import find_executable
from io import StringIO

import openmm
import openmm.app as app
import openmm.unit as u
import parmed
import yaml
from lxml import etree as et
from parmed.exceptions import ParameterWarning

warnings.filterwarnings("error", category=ParameterWarning)

_loadoffre = re.compile(r"loadoff (\S*)", re.I)
_sourcere = re.compile(r"source (\S*)", re.I)

# check for AMBERHOME, find from tleap location if not set, exception if can't
if os.getenv("AMBERHOME"):
    AMBERHOME = os.getenv("AMBERHOME")
else:
    if not find_executable("tleap"):
        raise Exception("AMBERHOME not set and tleap not available from PATH")
    tleap_path = find_executable("tleap")
    AMBERHOME = os.path.split(tleap_path)[0]
    AMBERHOME = os.path.join(AMBERHOME, "../")
    parmed.amber.AMBERHOME = AMBERHOME

# set global defaults for verbose and log
verbose = False
no_log = False

# set files that are ignored in leaprc's
# solvents and ions converted separately; leaprc.ff10 calls phosphoaa10.lib
# which does not exist anymore, LeAP skips it on error so we do too
ignore = {
    "solvents.lib",
    "atomic_ions.lib",
    "ions94.lib",
    "ions91.lib",
    "phosphoaa10.lib",
}

# define NEARLYZERO to replace numerical comparisons to zero
NEARLYZERO = 1e-10

# set beta
temperature = 300.0 * u.kelvin
kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
kT = kB * temperature
beta = 1.0 / kT


class LeapException(Exception):
    def __init__(self, leaprc_filename):
        msg = "Something went wrong in processing this LEaP input file:\n"
        msg += "\n"
        infile = open(leaprc_filename)
        contents = infile.read()
        msg += contents
        msg += "\n"
        super().__init__(msg)


def main():
    global verbose
    global no_log
    global logger
    # argparse
    parser = argparse.ArgumentParser(description="AMBER --> OpenMM forcefield " "conversion script")
    parser.add_argument(
        "--input",
        "-i",
        default="master.yaml",
        help='path of the input file. Default: "master.yaml"',
    )
    parser.add_argument(
        "--input-format",
        "-if",
        default="yaml",
        help='format of the input file: "yaml" or "leaprc". Default: "yaml"',
    )
    parser.add_argument(
        "--output-dir",
        "-od",
        help="path of the output directory. " 'Default: "ffxml/" for yaml, "./" for leaprc',
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="turns verbosity on")
    parser.add_argument(
        "--log",
        action="store",
        dest="log_filename",
        default=None,
        help="log energies for tests to specified CSV file",
    )
    parser.add_argument(
        "--protein-test",
        action="store_true",
        help="validate resulting XML through protein tests",
    )
    parser.add_argument(
        "--nucleic-test",
        action="store_true",
        help="validate resulting XML through nucleic acid tests",
    )
    parser.add_argument(
        "--protein-ua-test",
        action="store_true",
        help="validate resulting XML through united-atom protein tests",
    )
    parser.add_argument(
        "--phospho-protein-test",
        action="store_true",
        help="validate resulting XML through phosphorylated protein tests",
    )
    parser.add_argument(
        "--gaff-test",
        action="store_true",
        help="validate resulting XML through small-molecule (GAFF) test",
    )
    parser.add_argument(
        "--lipids-test",
        action="store_true",
        help="validate resulting XML through lipids tests",
    )
    args = parser.parse_args()
    verbose = args.verbose

    if args.log_filename:
        logger = Logger(args.log_filename)  # log to file
    else:
        logger = Logger()  # be silent

    # input is either a YAML or a leaprc - default is leaprc
    # output directory hardcoded here for ffxml/
    if args.input_format == "yaml":
        if args.output_dir is None:
            convert_yaml(args.input, ffxml_dir="../openmmforcefields/ffxml/amber")
        else:
            convert_yaml(args.input, ffxml_dir=args.output_dir)
    # if leaprc converted - output to the same dir
    elif args.input_format == "leaprc":
        if args.output_dir is None:
            ffxml_name = convert_leaprc(args.input, ffxml_dir="./")
        else:
            ffxml_name = convert_leaprc(args.input, ffxml_dir=args.output_dir)
        if args.protein_test:
            validate_protein(ffxml_name, args.input)
        if args.nucleic_test:
            validate_nucleic(ffxml_name, args.input)
        if args.protein_ua_test:
            validate_protein(ffxml_name, args.input, united_atom=True)
        if args.phospho_protein_test:
            validate_phospho_protein(ffxml_name, args.input)
        if args.gaff_test:
            validate_gaff(ffxml_name, args.input)
        if args.lipids_test:
            validate_lipids(ffxml_name, args.input)
    else:
        sys.exit("Wrong input_format chosen.")

    logger.close()


def read_lines(filename):
    """
    Read lines from file, stripping comments and newlines.

    """
    with open(filename) as f:
        lines = [line if "#" not in line else line[: line.index("#")] for line in f.readlines()]
    return lines


def write_file(file, contents):
    """
    Write text to file.

    Parameters
    ----------
    filename : str
       The file to write to
    contents : str
       Text contents to be written to file

    """
    if isinstance(file, str):
        outfile = open(file, "w")
    else:
        outfile = os.fdopen(file, "w")
    outfile.write(contents)
    outfile.close()


def convert_leaprc(
    files,
    split_filename=False,
    ffxml_dir="./",
    ignore=ignore,
    provenance=None,
    write_unused=False,
    override_level=None,
    filter_warnings="error",
    is_glycam=False,
):
    if verbose:
        print(f"\nConverting {files} to ffxml...")
    # allow for multiple source files - further code assuming list is passed
    if not isinstance(files, list):
        files = [files]
    basename = ""
    for f in files:
        f_basename = os.path.basename(f)
        if split_filename:
            f_basename = f_basename.split(".")[1:]
            f_basename = ".".join(f_basename)
        if not basename:
            basename = f_basename
        else:
            basename += "_"
            basename += f_basename
    ffxml_name = os.path.join(ffxml_dir, (basename + ".xml"))
    if not os.path.exists(ffxml_dir):
        os.mkdir(ffxml_dir)
    if verbose:
        print(f"Preprocessing the leaprc for {basename:s}...")
    # do source processing
    new_files = []
    for fil in files:
        lines = read_lines(fil)
        for line in lines:
            if _sourcere.findall(line):
                replace_leaprc = _sourcere.findall(line)[0]
                replace_leaprc_path = os.path.join(os.path.join(AMBERHOME, "dat/leap/cmd", replace_leaprc))
                new_files.append(replace_leaprc_path)
        new_files.append(fil)
    # now do ignore processing and join multiple leaprc's
    files = new_files
    new_lines = []
    for fil in files:
        lines = read_lines(fil)
        fil_new_lines = []
        for line in lines:
            if ignore is not None and _loadoffre.findall(line) and _loadoffre.findall(line)[0] in ignore:
                continue
            fil_new_lines += line
        new_lines += fil_new_lines
    leaprc = StringIO("".join(new_lines))
    if verbose:
        print(f"Converting to ffxml {ffxml_name}...")
    if filter_warnings != "error":
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
            params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, remediate_residues=(not write_unused))
    else:
        params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
        params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, remediate_residues=(not write_unused))
    if override_level:
        for name, residue in params.residues.items():
            residue.override_level = override_level
    if is_glycam:
        skip_duplicates = False
    else:
        skip_duplicates = True
    if filter_warnings != "error":
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(
                ffxml_name,
                provenance=provenance,
                write_unused=write_unused,
                improper_dihedrals_ordering="amber",
                skip_duplicates=skip_duplicates,
                is_glycam=is_glycam,
            )
    else:
        params.write(
            ffxml_name,
            provenance=provenance,
            write_unused=write_unused,
            improper_dihedrals_ordering="amber",
            skip_duplicates=skip_duplicates,
            is_glycam=is_glycam,
        )
    if verbose:
        print(f"{ffxml_name} successfully written!")
    return ffxml_name


def convert_gaff(
    files,
    ffxml_basename="",
    split_filename=False,
    ffxml_dir="../openmmforcefields/ffxml",
    ignore=ignore,
    provenance=None,
    write_unused=False,
    filter_warnings="error",
):
    if verbose:
        print(f"\nConverting {files} to ffxml...")
    # allow for multiple source files - further code assuming list is passed
    if not isinstance(files, list):
        files = [files]
    # Create ffxml
    ffxml_name = os.path.join(ffxml_dir, (ffxml_basename + ".xml"))
    if not os.path.exists(ffxml_dir):
        os.mkdir(ffxml_dir)
    # Process parameter file
    params = parmed.amber.AmberParameterSet(files)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, remediate_residues=(not write_unused))
    if filter_warnings != "error":
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(
                ffxml_name,
                provenance=provenance,
                write_unused=write_unused,
                improper_dihedrals_ordering="amber",
            )
    else:
        params.write(
            ffxml_name,
            provenance=provenance,
            write_unused=write_unused,
            improper_dihedrals_ordering="amber",
        )
    if verbose:
        print(f"{ffxml_name} successfully written!")
    return ffxml_name


def convert_recipe(
    files,
    solvent_file=None,
    ffxml_dir="./",
    provenance=None,
    ffxml_basename=None,
    filter_warnings="always",
):
    if verbose:
        print(f"\nConverting {files} to ffxml...")
    ffxml_name = os.path.join(ffxml_dir, (ffxml_basename + ".xml"))
    ffxml_temp_stringio = StringIO()
    params = parmed.amber.AmberParameterSet(files)
    print(params.atom_types.keys())
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    # Change atom type naming
    # atom_types
    new_atom_types = dict()
    for name, atom_type in params.atom_types.items():
        new_name = ffxml_basename + "-" + name
        new_atom_types[new_name] = atom_type
    params.atom_types = new_atom_types
    # atoms in residues
    for name, residue in params.residues.items():
        for atom in residue:
            new_type = ffxml_basename + "-" + atom.type
            atom.type = new_type
    if solvent_file is None:
        # this means this file does not include a water model - hard-coded assumption it is
        # then a 'multivalent' file - set overrideLevel to 1 for all residue templates
        for name, residue in params.residues.items():
            residue.override_level = 1
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(
                ffxml_name,
                provenance=provenance,
                write_unused=False,
                improper_dihedrals_ordering="amber",
            )
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(
                ffxml_temp_stringio,
                provenance=provenance,
                write_unused=False,
                improper_dihedrals_ordering="amber",
            )
        ffxml_temp_stringio.seek(0)
        if verbose:
            print("Modifying converted ffxml to append solvent parameters")
        tree_main = et.parse(ffxml_temp_stringio)
        tree_water = et.parse(solvent_file)
        root_main = tree_main.getroot()
        root_water = tree_water.getroot()
        with open(ffxml_name, "wb") as f:
            f.write(b"<ForceField>\n ")
            f.write(et.tostring(root_main.findall("Info")[0]))
            f.write(b"<AtomTypes>\n  ")
            for subelement in root_main.findall("AtomTypes")[0]:
                f.write(et.tostring(subelement))
            f.write(b" ")
            for subelement in root_water.findall("AtomTypes")[0]:
                f.write(et.tostring(subelement))
            f.write(b"</AtomTypes>\n <Residues>\n  ")
            for subelement in root_main.findall("Residues")[0]:
                f.write(et.tostring(subelement))
            f.write(b" ")
            for subelement in root_water.findall("Residues")[0]:
                f.write(et.tostring(subelement))
            f.write(b"</Residues>\n <HarmonicBondForce>\n  ")
            for subelement in root_water.findall("HarmonicBondForce")[0]:
                f.write(et.tostring(subelement))
            f.write(b"</HarmonicBondForce>\n <HarmonicAngleForce>\n  ")
            for subelement in root_water.findall("HarmonicAngleForce")[0]:
                f.write(et.tostring(subelement))
            f.write(b"</HarmonicAngleForce>\n ")
            f.write(
                (
                    '<NonbondedForce coulomb14scale="{}" lj14scale="{}">\n  '.format(
                        root_main.findall("NonbondedForce")[0].attrib["coulomb14scale"],
                        root_main.findall("NonbondedForce")[0].attrib["lj14scale"],
                    )
                ).encode("utf-8")
            )
            for subelement in root_main.findall("NonbondedForce")[0]:
                f.write(et.tostring(subelement))
            f.write(b" ")
            for subelement in root_water.findall("NonbondedForce")[0]:
                if subelement.tag == "UseAttributeFromResidue":
                    continue
                f.write(et.tostring(subelement))
            f.write(b"</NonbondedForce>\n</ForceField>")
    if verbose:
        print(f"{ffxml_name} successfully written!")
    return ffxml_name


def convert_yaml(yaml_name, ffxml_dir, ignore=ignore):
    data = yaml.load(open(yaml_name), Loader=yaml.FullLoader)
    # TODO: Verify that the version that is installed via conda matches sourcePackageVersion

    # Default yaml reading mode is leaprc
    ALLOWED_MODES = ("LEAPRC", "RECIPE", "GAFF")
    for entry in data:
        # Handle MODE switching
        if "MODE" in entry:
            MODE = entry["MODE"]
            if MODE not in ALLOWED_MODES:
                raise Exception(f"MODE definition must be one of {ALLOWED_MODES}")
            continue

        # Handle definition of source packages
        if "sourcePackage" in entry:
            source_pack = entry["sourcePackage"]
            source_pack_ver = entry["sourcePackageVersion"]
            continue
        if "sourcePackage2" in entry:
            # Switch mode to RECIPE processing
            source_pack2 = entry["sourcePackage2"]
            source_pack_ver2 = entry["sourcePackageVersion2"]
            continue

        # Extract source files, reference, and test files
        source_files = entry["Source"]
        reference = entry["Reference"]
        test_filename = entry["Test"]

        # Make sure source_files is a list
        if isinstance(source_files, str):
            source_files = [source_files]

        # Recipes require extra definitions
        if MODE == "RECIPE":
            recipe_name = entry["Name"]
            solvent_name = entry["Solvent"]
            if "Solvent_source" in entry:
                recipe_source2 = entry["Solvent_source"]
            else:
                recipe_source2 = None
            if "Standard" in entry:
                standard_ffxml = os.path.join(ffxml_dir, (entry["Standard"] + ".xml"))
            else:
                standard_ffxml = None
        elif MODE == "GAFF":
            recipe_name = entry["Name"]

        # Create provenance object
        provenance = dict()
        files = []
        source = provenance["Source"] = []
        for source_file in source_files:
            if MODE == "LEAPRC":
                if os.path.exists(source_file):
                    _filename = os.path.join("./", source_file)
                else:
                    _filename = os.path.join(AMBERHOME, "dat/leap/cmd", source_file)
            elif MODE == "RECIPE":
                _filename = os.path.join(AMBERHOME, "dat/leap/", source_file)
            elif MODE == "GAFF":
                _filename = os.path.join("../openmmforcefields/ffxml/amber/gaff/dat", source_file)
            files.append(_filename)
            source.append(dict())
            source[-1]["Source"] = source_file
            md5 = hashlib.md5()
            with open(_filename, "rb") as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]["md5hash"] = md5
            source[-1]["sourcePackage"] = source_pack
            source[-1]["sourcePackageVersion"] = source_pack_ver

        # For recipes, add water file and source info for it
        if MODE == "RECIPE" and recipe_source2 is not None:
            _filename = os.path.join("files", recipe_source2)
            solvent_file = _filename
            source.append(dict())
            source[-1]["Source"] = recipe_source2
            md5 = hashlib.md5()
            with open(_filename, "rb") as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]["md5hash"] = md5
            source[-1]["sourcePackage"] = source_pack2
            source[-1]["sourcePackageVersion"] = source_pack_ver2
        elif MODE == "RECIPE" and recipe_source2 is None:
            solvent_file = None
        provenance["Reference"] = reference

        # set default conversion options
        write_unused = False
        filter_warnings = "error"
        override_level = None
        # set conversion options if present
        if "Options" in entry:
            for option in entry["Options"]:
                if option == "write_unused":
                    write_unused = entry["Options"][option]
                elif option == "filter_warnings":
                    filter_warnings = entry["Options"][option]
                elif option == "ffxml_dir":
                    ffxml_dir = entry["Options"][option]
                elif option == "override_level":
                    override_level = entry["Options"][option]
                else:
                    raise Exception(f"Wrong option used in Options for {source_files:s}")

        # Convert files
        if MODE == "LEAPRC":
            is_glycam = False
            for source_file in source_files:
                if "GLYCAM" in source_file:
                    is_glycam = True
            ffxml_name = convert_leaprc(
                files,
                ffxml_dir=ffxml_dir,
                ignore=ignore,
                provenance=provenance,
                write_unused=write_unused,
                override_level=override_level,
                filter_warnings=filter_warnings,
                split_filename=True,
                is_glycam=is_glycam,
            )
        elif MODE == "RECIPE":
            ffxml_name = convert_recipe(
                files,
                solvent_file=solvent_file,
                ffxml_dir=ffxml_dir,
                provenance=provenance,
                ffxml_basename=recipe_name,
            )
        elif MODE == "GAFF":
            ffxml_name = convert_gaff(
                files,
                ffxml_basename=recipe_name,
                ffxml_dir=ffxml_dir,
                ignore=ignore,
                provenance=provenance,
                write_unused=write_unused,
                filter_warnings=filter_warnings,
                split_filename=True,
            )

        if "CharmmFFXMLFilename" in entry:
            charmm_ffxml_filename = entry["CharmmFFXMLFilename"]
            charmm_lipid2amber_filename = entry["CharmmLipid2AmberFilename"]
            merged_ffxml_filename = os.path.splitext(ffxml_name)[0] + "_merged.xml"
            if verbose:
                print("Merging lipid entries...")
            merge_lipids(ffxml_name, charmm_ffxml_filename, charmm_lipid2amber_filename, merged_ffxml_filename)
        if "Prefix" in entry:
            prefix = entry["Prefix"]
            if verbose:
                print(f'Rewriting {ffxml_name} to append prefix "{prefix}"...')
            add_prefix_to_ffxml(ffxml_name, prefix)
        if verbose:
            print("Validating the conversion...")
        tested = False
        for test in test_filename:
            if test == "protein":
                validate_protein(ffxml_name, entry["Source"])
                tested = True
            elif test == "nucleic":
                validate_dna(ffxml_name, entry["Source"])
                validate_rna(ffxml_name, entry["Source"])
                tested = True
            elif test == "protein_ua":
                validate_protein(ffxml_name, entry["Source"], united_atom=True)
                tested = True
            elif test == "protein_phospho":
                validate_phospho_protein(ffxml_name, entry["Source"])
                tested = True
            elif test == "gaff":
                validate_gaff(ffxml_name, entry["leaprc"], entry["Source"])
                tested = True
            elif test == "water_ion":
                validate_water_ion(
                    ffxml_name,
                    files,
                    solvent_name,
                    recipe_name,
                    standard_ffxml=standard_ffxml,
                )
                tested = True
            elif test == "dna":
                validate_dna(ffxml_name, entry["Source"])
                tested = True
            elif test == "rna":
                validate_rna(ffxml_name, entry["Source"])
                tested = True
            elif test == "lipids":
                validate_lipids(ffxml_name, entry["Source"])
                if "CharmmFFXMLFilename" in entry:
                    validate_merged_lipids(merged_ffxml_filename, entry["Source"])
                tested = True
            elif test == "protein_glycan":
                validate_glyco_protein(ffxml_name, entry["Source"])
                tested = True
        if not tested:
            raise Exception(f"No validation tests have been run for {source_files}")


def merge_lipids(ffxml_filename, charmm_ffxml_filename, charmm_lipid2amber_filename, merged_ffxml_filename):
    """
    Merge lipid residue definitions in AMBER ffxml file according to entries in a CHARMM ffxml file.

    Parameters
    ----------
    ffxml_filename : str
       AMBER lipids ffxml filename with AMBER split lipids.
    charmm_ffxml_filename : str
       CHARMM ffxml lipids
    charmmlipid2amber_filename : str
       CHARMM CSV file containing translation from CHARMM -> AMBER
    merged_ffxml_filename : str
       AMBER lipids ffxml filename merged with CHARMM lipids.
    """
    # Read the input files.
    charmmff = etree.parse(charmm_ffxml_filename)
    amberff = etree.parse(ffxml_filename)
    charmmResidues = charmmff.getroot().find("Residues").findall("Residue")
    amberResidues = amberff.getroot().find("Residues").findall("Residue")
    amberResMap = {}
    for res in amberResidues:
        atoms = {atom.attrib["name"]: atom for atom in res.findall("Atom")}
        amberResMap[res.attrib["name"]] = atoms
    translations = {}
    with open(charmm_lipid2amber_filename) as input:
        # Skip the first two lines.
        input.readline()
        input.readline()
        for line in input:
            fields = line.split(",")
            mergedRes = fields[0]
            mergedAtom = fields[2].split()[0]
            originalAtom, originalRes = fields[3].split()
            translations[(mergedRes, mergedAtom)] = (originalRes, originalAtom)

    # Remove all residues from the Amber file.
    parentNode = amberff.getroot().find("Residues")
    for res in amberResidues:
        parentNode.remove(res)

    # Copy over the CHARMM residues, making appropriate replacements.
    def translateResidue(residue):
        newres = deepcopy(residue)
        # Translate atom properties
        for atom in newres.findall("Atom"):
            key = (residue.attrib["name"], atom.attrib["name"])
            if key not in translations:
                return None  # We don't have a translation.
            amberResName, amberAtomName = translations[key]
            if amberResName not in amberResMap or amberAtomName not in amberResMap[amberResName]:
                return None  # We don't have a translation.
            amberAtom = amberResMap[amberResName][amberAtomName]
            for attrib in amberAtom.attrib:
                if attrib != "name":
                    atom.attrib[attrib] = amberAtom.attrib[attrib]
        # Remove Patches from CHARMM residues
        for patch in newres.findall("AllowPatch"):
            newres.remove(patch)

        return newres

    # Iterate over CHARMM lipid residue templates and replace components with AMBER parameters
    for residue in charmmResidues:
        copy = translateResidue(residue)
        if copy is not None:
            parentNode.append(copy)

    # Write merged lipid ffxml file (overwriting original file)
    amberff.write(merged_ffxml_filename)


def add_prefix_to_ffxml(ffxml_filename, prefix):
    """
    Replace the contents of an ffxml file with a modified version in which every atom type is prefixed with `prefix`.

    Parameters
    ----------
    ffxml_filename : str
       OpenMM ffxml filename (will be overwritten)
    prefix : str
       Prefix

    """

    import re

    inTypes = False
    replacements = {}

    modified_contents = ""
    with open(ffxml_filename) as infile:
        for line in infile:
            if "<AtomTypes>" in line:
                inTypes = True
            if "</AtomTypes>" in line:
                inTypes = False
            if inTypes:
                match = re.search('name="(.*?)"', line)
                if match is not None:
                    name = match.group(1)
                    newName = prefix + "-" + name
                    line = line.replace(f'name="{name}"', f'name="{newName}"')
                    replacements[f'type="{name:s}"'] = f'type="{newName}"'
                    replacements[f'type1="{name:s}"'] = f'type1="{newName}"'
                    replacements[f'type2="{name:s}"'] = f'type2="{newName}"'
                    replacements[f'type3="{name:s}"'] = f'type3="{newName}"'
                    replacements[f'type4="{name:s}"'] = f'type4="{newName}"'
            else:
                for key in replacements:
                    if key in line:
                        line = line.replace(key, replacements[key])
            if line.endswith("\n"):
                line = line[:-1]
            modified_contents += line + "\n"

    with open(ffxml_filename, "w") as outfile:
        outfile.write(modified_contents)


def assert_energies_glyco_protein(prmtop, inpcrd, ffxml, tolerance=1e-1):
    import math

    # Get AMBER system
    parm_amber = parmed.load_file(prmtop, inpcrd)
    system_amber = parm_amber.createSystem()

    # Create topology where residue names are named from HYP to CHYP or NHYP (etc) where necessary
    source_topology = parm_amber.topology
    destination_topology = app.Topology()

    new_atoms = {}
    for chain in source_topology.chains():
        new_chain = destination_topology.addChain(chain.id)
        for residue in chain.residues():
            new_name = residue.name
            if residue.index in [0, 5, 13, 21, 29]:
                new_name = "N" + residue.name
            elif residue.index in [4, 9, 17, 25, 33]:
                new_name = "C" + residue.name
            new_residue = destination_topology.addResidue(new_name, new_chain, residue.id)
            for atom in residue.atoms():
                new_atom = destination_topology.addAtom(atom.name, atom.element, new_residue, atom.id)
                new_atoms[atom] = new_atom
    for bond in source_topology.bonds():
        order = bond.order
        destination_topology.addBond(new_atoms[bond[0]], new_atoms[bond[1]], order=order)

    # Get OpenMM system
    if isinstance(ffxml, str):
        ff = app.ForceField(ffxml)
    else:
        ff = app.ForceField(*ffxml)
    system_omm = ff.createSystem(destination_topology, ignoreExternalBonds=True)

    def compute_potential_components(system, positions, beta=beta):
        # Note: this is copied from perses
        system = deepcopy(system)
        for index in range(system.getNumForces()):
            force = system.getForce(index)
            force.setForceGroup(index)
        integrator = openmm.VerletIntegrator(1.0 * u.femtosecond)
        platform = openmm.Platform.getPlatformByName("Reference")
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        energy_components = dict()
        for index in range(system.getNumForces()):
            force = system.getForce(index)
            forcename = force.__class__.__name__
            groups = 1 << index
            potential = beta * context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
            energy_components.setdefault(forcename, 0.0)
            energy_components[forcename] += potential
        del context, integrator
        return energy_components

    amber_energies = compute_potential_components(system_amber, parm_amber.positions)
    omm_energies = compute_potential_components(system_omm, parm_amber.positions)

    for force_name in sorted(set(amber_energies) | set(omm_energies)):
        assert math.isclose(amber_energies.get(force_name, 0.0), omm_energies.get(force_name, 0.0), rel_tol=tolerance)


def assert_energies(
    prmtop,
    inpcrd,
    ffxml,
    system_name="unknown",
    tolerance=2.5e-5,
    improper_tolerance=1e-2,
    units=u.kilojoules_per_mole,
    openmm_topology=None,
    openmm_positions=None,
):
    # AMBER
    parm_amber = parmed.load_file(prmtop, inpcrd)
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=units)

    # OpenMM-ffxml
    if isinstance(ffxml, str):
        ff = app.ForceField(ffxml)
    else:
        ff = app.ForceField(*ffxml)

    if openmm_positions is None:
        openmm_positions = parm_amber.positions

    if openmm_topology is not None:
        system_omm = ff.createSystem(openmm_topology)
        parm_omm = parmed.openmm.load_topology(openmm_topology, system_omm, xyz=openmm_positions)
    else:
        system_omm = ff.createSystem(parm_amber.topology)
        parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    
    try:
        system_omm = parm_omm.createSystem(splitDihedrals=True)
    except parmed.exceptions.ParameterError:
        # Unfortunately ParmEd does not understand how to interpret the systems
        # it creates with NBFIX/LJEDIT.  To handle this, first get the nonbonded
        # energies from the two systems.
        nonbonded_names = ("NonbondedForce", "CustomNonbondedForce", "CustomBondForce")
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=units, platform="Reference")
        amber_energy_nb = sum(energy for name, energy in amber_energies if name in nonbonded_names)
        omm_energy_nb = sum(energy for name, energy in amber_energies if name in nonbonded_names)
        
        # Delete all of the nonbonded forces from the OpenMM system (that has
        # not been sent through ParmEd yet).
        for force_index in range(system_omm.getNumForces() - 1, -1, -1):
            if isinstance(system_omm.getForce(force_index), (openmm.NonbondedForce, openmm.CustomBondForce, openmm.CustomNonbondedForce)):
                system_omm.removeForce(force_index)

        # Send the system round-trip through ParmEd to split the dihedrals and
        # re-evaluate the energies.
        parm_omm = parmed.openmm.load_topology(parm_amber.topology if openmm_topology is None else openmm_topology, system_omm, xyz=openmm_positions)
        system_omm = parm_omm.createSystem(splitDihedrals=True)
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=units, platform="Reference")

        # Remove the nonbonded energies from both lists, and add back the ones
        # we manually calculated above.
        amber_energies = [(name, energy) for name, energy in amber_energies if name not in nonbonded_names]
        omm_energies = [(name, energy) for name, energy in omm_energies if name not in nonbonded_names]
        amber_energies.append(("NonbondedForce", amber_energy_nb))
        omm_energies.append(("NonbondedForce", omm_energy_nb))

    else:
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=units, platform="Reference")

    # calc rel energies and assert
    _energies = []
    rel_energies = []
    for i, j in zip(amber_energies, omm_energies):
        if i[0] != j[0]:
            raise Exception("Mismatch in energy tuples naming.")
        if abs(i[1]) > NEARLYZERO:
            rel_energies.append((i[0], abs((i[1] - j[1]) / i[1])))
        else:
            if abs(j[1]) > NEARLYZERO:
                raise AssertionError(
                    f"One of AMBER {system_name} energies ({i[0]}) for {ffxml} is zero, "
                    "while the corresponding OpenMM energy is non-zero"
                )
            rel_energies.append((i[0], 0))

    dihedrals_done = False
    for i, amber_energy, openmm_energy in zip(rel_energies, amber_energies, omm_energies):
        if i[0] != "PeriodicTorsionForce":
            if i[1] > tolerance:
                raise AssertionError(
                    f"{system_name} relative energy error ({i[0]}, {i[1]:f}) "
                    f"outside of allowed tolerance ({tolerance:f}) for {ffxml}: "
                    f"AMBER {amber_energy} OpenMM {openmm_energy}"
                )
        else:
            if not dihedrals_done:
                if i[1] > tolerance:
                    raise AssertionError(
                        f"{system_name} relative energy error ({i[0]}, {i[1]:f}) "
                        f"outside of allowed tolerance ({tolerance:f}) for {ffxml}: "
                        f"AMBER {amber_energy} OpenMM {openmm_energy}"
                    )
                dihedrals_done = True
            else:  # impropers
                if i[1] > improper_tolerance:
                    raise AssertionError(
                        f"{system_name} relative energy error ({i[0]}-impropers, {i[1]:f}) "
                        f"outside of allowed tolerance ({improper_tolerance:f}) for {ffxml}: "
                        f"AMBER {amber_energy} OpenMM {openmm_energy}"
                    )

    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log["ffxml_name"] = ffxml
        amber_energies_log["test_system"] = system_name
        amber_energies_log["data_type"] = "AMBER"
        amber_energies_log["units"] = units
        omm_energies_log["ffxml_name"] = ffxml
        omm_energies_log["test_system"] = system_name
        omm_energies_log["data_type"] = "OpenMM"
        omm_energies_log["units"] = units
        rel_energies_log["ffxml_name"] = ffxml
        rel_energies_log["test_system"] = system_name
        rel_energies_log["data_type"] = "abs((AMBER-OpenMM)/AMBER)"
        dihedrals_done = False
        for item in amber_energies:
            if item[0] == "PeriodicTorsionForce" and not dihedrals_done:
                amber_energies_log["PeriodicTorsionForce_dihedrals"] = item[1]
                dihedrals_done = True
            elif item[0] == "PeriodicTorsionForce" and dihedrals_done:
                amber_energies_log["PeriodicTorsionForce_impropers"] = item[1]
            elif item[0] == "CMMotionRemover":
                continue
            else:
                amber_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in omm_energies:
            if item[0] == "PeriodicTorsionForce" and not dihedrals_done:
                omm_energies_log["PeriodicTorsionForce_dihedrals"] = item[1]
                dihedrals_done = True
            elif item[0] == "PeriodicTorsionForce" and dihedrals_done:
                omm_energies_log["PeriodicTorsionForce_impropers"] = item[1]
            elif item[0] == "CMMotionRemover":
                continue
            else:
                omm_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in rel_energies:
            if item[0] == "PeriodicTorsionForce" and not dihedrals_done:
                rel_energies_log["PeriodicTorsionForce_dihedrals"] = item[1]
                dihedrals_done = True
            elif item[0] == "PeriodicTorsionForce" and dihedrals_done:
                rel_energies_log["PeriodicTorsionForce_impropers"] = item[1]
            elif item[0] == "CMMotionRemover":
                continue
            else:
                rel_energies_log[item[0]] = item[1]

        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)


def validate_protein(ffxml_name, leaprc_name, united_atom=False):
    if verbose:
        print(f"Protein energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    ala3_top = tempfile.mkstemp()
    ala3_crd = tempfile.mkstemp()
    villin_top = tempfile.mkstemp()
    villin_crd = tempfile.mkstemp()
    leap_script_ala3_file = tempfile.mkstemp()
    leap_script_villin_file = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")
    if not united_atom:
        leap_script_ala3_string = f"""source {leaprc_name}
x = loadPdb files/ala3.pdb
saveAmberParm x {ala3_top[1]} {ala3_crd[1]}
quit"""
        leap_script_villin_string = f"""source {leaprc_name}
x = loadPdb files/villin.pdb
saveAmberParm x {villin_top[1]} {villin_crd[1]}
quit"""
    else:
        leap_script_ala3_string = f"""source {leaprc_name}
x = loadPdb files/ala3_ua.pdb
saveAmberParm x {ala3_top[1]} {ala3_crd[1]}
quit"""
        leap_script_villin_string = f"""source {leaprc_name}
x = loadPdb files/villin_ua.pdb
saveAmberParm x {villin_top[1]} {villin_crd[1]}
quit"""

    write_file(leap_script_ala3_file[0], leap_script_ala3_string)
    write_file(leap_script_villin_file[0], leap_script_villin_string)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_ala3_file[1]} > {os.devnull}")
    if os.path.getsize(ala3_top[1]) == 0 or os.path.getsize(ala3_crd[1]) == 0:
        raise LeapException(leap_script_ala3_file[1])
    os.system(f"tleap -f {leap_script_villin_file[1]} > {os.devnull}")
    if os.path.getsize(villin_top[1]) == 0 or os.path.getsize(villin_crd[1]) == 0:
        raise LeapException(leap_script_villin_file[1])

    try:
        if verbose:
            print("Calculating and validating ala_ala_ala energies...")
        assert_energies(ala3_top[1], ala3_crd[1], ffxml_name, system_name="protein-ala_ala_ala")
        if verbose:
            print("Ala_ala_ala energy validation successful!")

        if verbose:
            print("Calculating and validating villin headpiece energies...")
        assert_energies(
            villin_top[1],
            villin_crd[1],
            ffxml_name,
            system_name="protein-villin headpiece",
        )
        if verbose:
            print("Villin headpiece energy validation successful!")
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (
            ala3_top,
            ala3_crd,
            villin_top,
            villin_crd,
            leap_script_ala3_file,
            leap_script_villin_file,
        ):
            os.unlink(f[1])
    if verbose:
        print(f"Protein energy validation for {ffxml_name} done!")


def validate_dna(ffxml_name, leaprc_name):
    if verbose:
        print(f"DNA energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    dna_top = tempfile.mkstemp()
    dna_crd = tempfile.mkstemp()
    leap_script_dna_file = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_dna_string = f"""addPdbAtomMap {{
{{ "H1'" "H1*" }}
{{ "H2'" "H2'1" }}
{{ "H2''" "H2'2" }}
{{ "H3'" "H3*" }}
{{ "H4'" "H4*" }}
{{ "H5'" "H5'1" }}
{{ "H5''" "H5'2" }}
{{ "HO2'" "HO'2" }}
{{ "HO5'" "H5T"  }}
{{ "HO3'" "H3T" }}
{{ "OP1" "O1P" }}
{{ "OP2" "O2P" }}
}}
source {leaprc_name}
addPdbResMap {{
{{ 0 "DG" "DG5"  }} {{ 1 "DG" "DG3"  }}
{{ 0 "DA" "DA5"  }} {{ 1 "DA" "DA3"  }}
{{ 0 "DC" "DC5"  }} {{ 1 "DC" "DC3"  }}
{{ 0 "DT" "DT5"  }} {{ 1 "DT" "DT3"  }}
}}
x = loadPdb files/4rzn_dna.pdb
saveAmberParm x {dna_top[1]} {dna_crd[1]}
quit"""

    write_file(leap_script_dna_file[0], leap_script_dna_string)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_dna_file[1]} > {os.devnull}")
    if os.path.getsize(dna_top[1]) == 0 or os.path.getsize(dna_crd[1]) == 0:
        raise LeapException(leap_script_dna_file[1])

    try:
        if verbose:
            print("Calculating and validating DNA energies...")
        assert_energies(dna_top[1], dna_crd[1], ffxml_name, system_name="nucleic-DNA")
        if verbose:
            print("DNA energy validation successful!")

    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (dna_top, dna_crd, leap_script_dna_file):
            os.unlink(f[1])
    if verbose:
        print(f"DNA energy validation for {ffxml_name} done!")


def validate_rna(ffxml_name, leaprc_name):
    if verbose:
        print(f"RNA energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    rna_top = tempfile.mkstemp()
    rna_crd = tempfile.mkstemp()
    leap_script_rna_file = tempfile.mkstemp()
    leap_script_rna_file_alt = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")

    leap_script_rna_string = f"""
addPdbAtomMap {{
{{ "H1'" "H1*" }}
{{ "H2'" "H2'1" }}
{{ "H2''" "H2'2" }}
{{ "H3'" "H3*" }}
{{ "H4'" "H4*" }}
{{ "H5'" "H5'1" }}
{{ "H5''" "H5'2" }}
{{ "HO2'" "HO'2" }}
{{ "HO5'" "H5T"  }}
{{ "HO3'" "H3T" }}
{{ "OP1" "O1P" }}
{{ "OP2" "O2P" }}
}}
source {leaprc_name}
addPdbResMap {{
{{ 0 "G" "G5"  }} {{ 1 "G" "G3"  }} {{ "G" "G" }}
{{ 0 "A" "A5"  }} {{ 1 "A" "A3"  }} {{ "A" "A" }}
{{ 0 "C" "C5"  }} {{ 1 "C" "C3"  }} {{ "C" "C" }}
{{ 0 "U" "U5"  }} {{ 1 "U" "U3"  }} {{ "U" "U" }}
}}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x {rna_top[1]} {rna_crd[1]}
quit"""

    leap_script_rna_string_alt = f"""
addPdbAtomMap {{
{{ "H1'" "H1*" }}
{{ "H2'" "H2'1" }}
{{ "H2''" "H2'2" }}
{{ "H3'" "H3*" }}
{{ "H4'" "H4*" }}
{{ "H5'" "H5'1" }}
{{ "H5''" "H5'2" }}
{{ "HO2'" "HO'2" }}
{{ "HO5'" "H5T"  }}
{{ "HO3'" "H3T" }}
{{ "OP1" "O1P" }}
{{ "OP2" "O2P" }}
}}
source {leaprc_name}
addPdbResMap {{
{{ 0 "G" "RG5"  }} {{ 1 "G" "RG3"  }} {{ "G" "RG" }}
{{ 0 "A" "RA5"  }} {{ 1 "A" "RA3"  }} {{ "A" "RA" }}
{{ 0 "C" "RC5"  }} {{ 1 "C" "RC3"  }} {{ "C" "RC" }}
{{ 0 "U" "RU5"  }} {{ 1 "U" "RU3"  }} {{ "U" "RU" }}
}}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x {rna_top[1]} {rna_crd[1]}
quit"""

    write_file(leap_script_rna_file[0], leap_script_rna_string)
    write_file(leap_script_rna_file_alt[0], leap_script_rna_string_alt)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_rna_file[1]} > {os.devnull}")
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        # try alternative name mappings
        os.system(f"tleap -f {leap_script_rna_file_alt[1]} > {os.devnull}")
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        raise LeapException(leap_script_rna_file_alt[1])

    try:
        if verbose:
            print("Calculating and validating RNA energies...")
        # improper testing turned off pending solution to problems
        assert_energies(rna_top[1], rna_crd[1], ffxml_name, system_name="nucleic-RNA")
        if verbose:
            print("RNA energy validation successful!")
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (rna_top, rna_crd, leap_script_rna_file, leap_script_rna_file_alt):
            os.unlink(f[1])
    if verbose:
        print(f"RNA energy validation for {ffxml_name} done!")


def validate_gaff(ffxml_name, leaprc_name, gaff_dat_name):
    if verbose:
        print(f"GAFF energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    imatinib_top = tempfile.mkstemp()
    imatinib_crd = tempfile.mkstemp()
    leap_script_imatinib_file = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_imatinib_string = f"""\
source {leaprc_name}
loadamberparams ../openmmforcefields/ffxml/amber/gaff/dat/{gaff_dat_name}
loadamberparams files/frcmod.imatinib
x = loadMol2 files/imatinib.mol2
saveAmberParm x {imatinib_top[1]} {imatinib_crd[1]}
quit"""
    write_file(leap_script_imatinib_file[0], leap_script_imatinib_string)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_imatinib_file[1]} > {os.devnull}")
    if os.path.getsize(imatinib_top[1]) == 0 or os.path.getsize(imatinib_crd[1]) == 0:
        raise LeapException(leap_script_imatinib_file[1])

    try:
        if verbose:
            print("Calculating and validating imatinib energies...")
        assert_energies(
            imatinib_top[1],
            imatinib_crd[1],
            (ffxml_name, "files/imatinib.xml", "files/imatinib_frcmod.xml"),
            system_name="gaff-imatinib",
        )
        if verbose:
            print("Imatinib energy validation successful!")
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (imatinib_top, imatinib_crd, leap_script_imatinib_file):
            os.unlink(f[1])
    if verbose:
        print(f"GAFF energy validation for {ffxml_name} done!")


def validate_phospho_protein(
    ffxml_name,
    leaprc_name,
    supp_leaprc_name="oldff/leaprc.ff99SBildn",
    supp_ffxml_name="../openmmforcefields/ffxml/amber/ff99SBildn.xml",
    phospho="phospho10",
):
    if "14" in leaprc_name:
        # Use AMBER14SB
        supp_leaprc_name = "oldff/leaprc.ff14SB"
        supp_ffxml_name = "../openmmforcefields/ffxml/amber/ff14SB.xml"
        phospho = "phospho14"

    # this function assumes ffxml/amber/ff14SB.xml already exists
    if verbose:
        print(f"Phosphorylated protein energy validation for {ffxml_name}")
    for pdbname in glob.iglob(f"files/{phospho}/*.pdb"):
        if verbose:
            print(f"Now testing with pdb {os.path.basename(pdbname)}")
        if verbose:
            print("Preparing temporary files for validation...")
        top = tempfile.mkstemp()
        crd = tempfile.mkstemp()
        leap_script_file = tempfile.mkstemp()

        if verbose:
            print("Preparing LeaP scripts...")
        leap_script_string = f"""source {supp_leaprc_name}
source {leaprc_name}
x = loadPdb {pdbname}
saveAmberParm x {top[1]} {crd[1]}
quit"""

        write_file(leap_script_file[0], leap_script_string)

        if verbose:
            print("Running LEaP...")
        os.system(f"tleap -f {leap_script_file[1]} > {os.devnull}")
        if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
            raise LeapException(leap_script_file[1])

        try:
            if verbose:
                print("Calculating and validating energies...")
            assert_energies(
                top[1],
                crd[1],
                (supp_ffxml_name, ffxml_name),
                system_name=f"phospho_protein: {os.path.basename(pdbname)}",
            )
            if verbose:
                print("Energy validation successful!")
        finally:
            if verbose:
                print("Deleting temp files...")
            for f in (top, crd, leap_script_file):
                os.unlink(f[1])
        if verbose:
            print(f"Phosphorylated protein energy validation for {ffxml_name} done!")


def validate_water_ion(ffxml_name, source_recipe_files, solvent_name, recipe_name, standard_ffxml=None):
    if verbose:
        print(f"Water and ions energy validation for {ffxml_name}")
    if solvent_name == "tip3p":
        HOH = "TP3"
        solvent_frcmod = "frcmod.tip3p"
    elif solvent_name == "tip4pew":
        HOH = "T4E"
        solvent_frcmod = "frcmod.tip4pew"
    elif solvent_name == "spce":
        HOH = "SPC"
        solvent_frcmod = "frcmod.spce"
    elif solvent_name == "tip3pfb":
        HOH = "FB3"
        solvent_frcmod = "frcmod.tip3pfb"
    elif solvent_name == "tip4pfb":
        HOH = "FB4"
        solvent_frcmod = "frcmod.tip4pfb"
    pdb_name = "files/water_ion/" + recipe_name + ".pdb"
    if verbose:
        print("Preparing temporary files for validation...")
    top = tempfile.mkstemp()
    crd = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()
    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_string_part1 = f"""loadamberparams parm10.dat
loadamberparams {source_recipe_files[0]}
loadamberparams {source_recipe_files[1]}\n"""

    leap_script_string_part2 = f"""\nloadOff atomic_ions.lib
loadoff solvents.lib
HOH = {HOH}
# for TIP4PEW
addPdbAtomMap {{{{ "M" "EPW" }}}}
x = loadPdb {pdb_name}
saveAmberParm x {top[1]} {crd[1]}
quit"""

    if solvent_frcmod:
        leap_script_string = (
            leap_script_string_part1 + (f"loadamberparams {solvent_frcmod}") + leap_script_string_part2
        )
    else:
        leap_script_string = leap_script_string_part1 + leap_script_string_part2

    write_file(leap_script_file[0], leap_script_string)

    # this test does it's own energy assertion because of differences
    if verbose:
        print("Running LEaP...")

    os.system(f"tleap -f {leap_script_file[1]} > {os.devnull}")
    if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
        raise LeapException(leap_script_file[1])
    try:
        if verbose:
            print("Calculating and validating energies...")
        pdb = app.PDBFile(pdb_name, extraParticleIdentifier="")
        if standard_ffxml is None:
            ff = app.ForceField(ffxml_name)
        else:
            ff = app.ForceField(ffxml_name, standard_ffxml)
        system_omm = ff.createSystem(pdb.topology)
        parm_omm = parmed.openmm.load_topology(pdb.topology, xyz=pdb.positions)
        parm_amber = parmed.load_file(top[1], crd[1])
        system_amber = parm_amber.createSystem()
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)
        for entry in omm_energies:
            if entry[0] == "NonbondedForce":
                omm_nonbonded = entry[1]
        amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
        for entry in amber_energies:
            if entry[0] == "NonbondedForce":
                amber_nonbonded = entry[1]

        rel_nonbonded = abs((amber_nonbonded - omm_nonbonded) / amber_nonbonded)

        if rel_nonbonded > 1e-5:
            raise AssertionError(
                f"NonbondedForce Water and ions energy ({rel_nonbonded:f}) outside of "
                f"allowed tolerance ({1e-5:f}) for {ffxml_name}:"
            )
        if verbose:
            print("Energy validation successful!")

    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (top, crd, leap_script_file):
            os.unlink(f[1])
    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log["ffxml_name"] = ffxml_name
        amber_energies_log["test_system"] = "water_ion"
        amber_energies_log["data_type"] = "AMBER"
        amber_energies_log["NonbondedForce"] = amber_nonbonded
        amber_energies_log["units"] = u.kilojoules_per_mole
        omm_energies_log["ffxml_name"] = ffxml_name
        omm_energies_log["test_system"] = "water_ion"
        omm_energies_log["data_type"] = "OpenMM"
        omm_energies_log["NonbondedForce"] = omm_nonbonded
        omm_energies_log["units"] = u.kilojoules_per_mole
        rel_energies_log["ffxml_name"] = ffxml_name
        rel_energies_log["test_system"] = "water_ion"
        rel_energies_log["data_type"] = "abs((AMBER-OpenMM)/AMBER)"
        rel_energies_log["NonbondedForce"] = rel_nonbonded
        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)
    if verbose:
        print(f"Water and ions energy validation for {ffxml_name} done!")


def validate_impropers(ffxml_name, leaprc_name):
    if verbose:
        print(f"Impropers validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    top_villin = tempfile.mkstemp()
    crd_villin = tempfile.mkstemp()
    top_dna = tempfile.mkstemp()
    crd_dna = tempfile.mkstemp()
    top_rna = tempfile.mkstemp()
    crd_rna = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_string = f"""source {leaprc_name}
x = loadPdb files/villin.pdb
y = loadPdb files/4rzn_dna.pdb
z = loadPdb files/5c5w_rna.pdb
saveAmberParm x {top_villin[1]} {crd_villin[1]}
saveAmberParm y {top_dna[1]} {crd_dna[1]}
saveAmberParm z {top_rna[1]} {crd_rna[1]}
quit"""
    write_file(leap_script_file[0], leap_script_string)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_file[1]} > {os.devnull}")
    if os.path.getsize(top_villin[1]) == 0 or os.path.getsize(crd_villin[1]) == 0:
        raise LeapException(leap_script_file[1])
    if os.path.getsize(top_dna[1]) == 0 or os.path.getsize(crd_dna[1]) == 0:
        raise LeapException(leap_script_file[1])
    if os.path.getsize(top_rna[1]) == 0 or os.path.getsize(crd_rna[1]) == 0:
        raise LeapException(leap_script_file[1])

    # load into parmed
    parm_amber_villin = parmed.load_file(top_villin[1])
    parm_amber_dna = parmed.load_file(top_dna[1])
    parm_amber_rna = parmed.load_file(top_rna[1])

    # OpenMM
    ff = app.ForceField(ffxml_name)
    sys_omm_villin = ff.createSystem(parm_amber_villin.topology)
    sys_omm_dna = ff.createSystem(parm_amber_dna.topology)
    sys_omm_rna = ff.createSystem(parm_amber_rna.topology)
    parm_omm_villin = parmed.openmm.load_topology(parm_amber_villin.topology, sys_omm_villin)
    parm_omm_dna = parmed.openmm.load_topology(parm_amber_dna.topology, sys_omm_dna)
    parm_omm_rna = parmed.openmm.load_topology(parm_amber_rna.topology, sys_omm_rna)

    # prepare sets of idxs
    set_amber_villin = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx)
        for dih in parm_amber_villin.dihedrals
        if dih.improper
    }
    set_amber_dna = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx)
        for dih in parm_amber_dna.dihedrals
        if dih.improper
    }
    set_amber_rna = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx)
        for dih in parm_amber_rna.dihedrals
        if dih.improper
    }
    set_omm_villin = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx)
        for dih in parm_omm_villin.dihedrals
        if dih.improper
    }
    set_omm_dna = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx) for dih in parm_omm_dna.dihedrals if dih.improper
    }
    set_omm_rna = {
        (dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx) for dih in parm_omm_rna.dihedrals if dih.improper
    }

    try:
        if set_amber_villin - set_omm_villin != set() or set_omm_villin - set_amber_villin != set():
            raise AssertionError(
                f"""Impropers validation fail for {ffxml_name} (villin)
                                    set_amber - set_omm: {set_amber_villin - set_omm_villin}
                                    set_omm - set_amber: {set_omm_villin - set_amber_villin}"""
            )
        if set_amber_dna - set_omm_dna != set() or set_omm_dna - set_amber_dna != set():
            raise AssertionError(
                f"""Impropers validation fail for {ffxml_name} (DNA)
                                    set_amber - set_omm: {set_amber_dna - set_omm_dna}
                                    set_omm - set_amber: {set_omm_dna - set_amber_dna}"""
            )
        if set_amber_rna - set_omm_rna != set() or set_omm_rna - set_amber_rna != set():
            raise AssertionError(
                f"""Impropers validation fail for {ffxml_name} (RNA)
                                    set_amber - set_omm: {set_amber_rna - set_omm_rna}
                                    set_omm - set_amber: {set_omm_rna - set_amber_rna}"""
            )
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (
            top_villin,
            crd_villin,
            top_dna,
            crd_dna,
            top_rna,
            crd_rna,
            leap_script_file,
        ):
            os.unlink(f[1])
    if verbose:
        print(f"Improper validation for {ffxml_name} done!")


def validate_lipids(ffxml_name, leaprc_name):
    if verbose:
        print(f"Lipids energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    lipids_top = tempfile.mkstemp()
    lipids_crd = tempfile.mkstemp()
    leap_script_lipids_file = tempfile.mkstemp()

    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_lipids_string = f"""source {leaprc_name}
x = loadPdb files/POPC-nowater-amber.pdb
saveAmberParm x {lipids_top[1]} {lipids_crd[1]}
quit"""
    write_file(leap_script_lipids_file[0], leap_script_lipids_string)

    if verbose:
        print("Running LEaP...")

    os.system(f"tleap -f {leap_script_lipids_file[1]} > {os.devnull}")
    if os.path.getsize(lipids_top[1]) == 0 or os.path.getsize(lipids_crd[1]) == 0:
        raise LeapException(leap_script_lipids_file[1])

    try:
        if verbose:
            print("Calculating and validating lipids energies...")
        assert_energies(lipids_top[1], lipids_crd[1], ffxml_name, system_name="lipids")
        if verbose:
            print("Lipids energy validation successful!")
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (lipids_top, lipids_crd, leap_script_lipids_file):
            os.unlink(f[1])
    if verbose:
        print(f"Lipids energy validation for {ffxml_name} done!")


def validate_merged_lipids(ffxml_name, leaprc_name):
    if verbose:
        print(f"Lipids (merged) energy validation for {ffxml_name}")
    if verbose:
        print("Preparing temporary files for validation...")
    lipids_top = tempfile.mkstemp()
    lipids_crd = tempfile.mkstemp()
    leap_script_lipids_file = tempfile.mkstemp()
    pdbfile = app.PDBFile("files/POPC-nowater-charmm.pdb")

    if verbose:
        print("Preparing LeaP scripts...")
    leap_script_lipids_string = f"""source {leaprc_name}
x = loadPdb files/POPC-nowater-amber.pdb
saveAmberParm x {lipids_top[1]} {lipids_crd[1]}
quit"""
    write_file(leap_script_lipids_file[0], leap_script_lipids_string)

    if verbose:
        print("Running LEaP...")
    os.system(f"tleap -f {leap_script_lipids_file[1]} > {os.devnull}")
    if os.path.getsize(lipids_top[1]) == 0 or os.path.getsize(lipids_crd[1]) == 0:
        raise LeapException(leap_script_lipids_file[1])

    try:
        if verbose:
            print("Calculating and validating lipids energies...")
        assert_energies(
            lipids_top[1],
            lipids_crd[1],
            ffxml_name,
            system_name="lipids",
            openmm_topology=pdbfile.topology,
            openmm_positions=pdbfile.positions,
        )
        if verbose:
            print("Lipids energy validation successful!")
    finally:
        if verbose:
            print("Deleting temp files...")
        for f in (lipids_top, lipids_crd, leap_script_lipids_file):
            os.unlink(f[1])
    if verbose:
        print(f"Lipids energy validation for {ffxml_name} done!")


def modify_glycan_ffxml(input_ffxml_path):
    """
    Creates a modified XML file with the following changes:
    - All parameters are specified by type, not class.
    - All atom type definitions that are not GLYCAM-specific are removed.
    - All parameter definitions that do not involve at least one GLYCAM-specific type are removed.
    - For standard (not GLYCAM-specific) atom types, it adds the protein- prefix.
    - Updates unscaled types in the script
    - The bad torsion (NH-Cg-Cg-Sm, bad because NH atom type is not defined) is removed.
    - GlycamTemplateMatcher (used for creating the system) initialization script is added.

    Makes the following assumptions about the input ffxml file:
    - No prefixes have been added.

    Parameters
    ----------
    input_ffxml_path : str
        path to ffxml for which to modify (in place)

    """

    import re
    import xml.etree.ElementTree as etree

    # Define atom types

    protein_types = {
        "C",
        "CA",
        "CB",
        "CC",
        "CN",
        "CR",
        "CT",
        "CV",
        "CW",
        "C*",
        "CX",
        "H",
        "HC",
        "H1",
        "HA",
        "H4",
        "H5",
        "HO",
        "HS",
        "HP",
        "N",
        "NA",
        "NB",
        "N2",
        "N3",
        "O",
        "O2",
        "OH",
        "S",
        "SH",
        "CO",
        "2C",
        "3C",
        "C8",
    }
    solvent_types = {
        "HW",
        "OW",
        "Li+",
        "Na+",
        "K+",
        "Rb+",
        "Cs+",
        "F-",
        "Cl-",
        "Br-",
        "I-",
    }
    glycam_types = set()
    replacements = {}
    for type in protein_types:
        replacements[type] = "protein-" + type

    # Process <AtomTypes>

    tree = etree.parse(input_ffxml_path)
    root = tree.getroot()
    types = root.find("AtomTypes")
    for type in types.findall("Type"):
        name = type.get("name")
        if name in protein_types or name in solvent_types:
            types.remove(type)
        else:
            glycam_types.add(name)
            replacements[name] = "glycam-" + name
            type.set("name", replacements[name])

    # Process <Residues>

    residues = root.find("Residues")
    for residue in residues.findall("Residue"):
        for atom in residue.findall("Atom"):
            atom.set("type", replacements[atom.get("type")])

    # Process <HarmonicBondForce>

    force = root.find("HarmonicBondForce")
    for bond in force.findall("Bond"):
        # Change attributes from class to type
        bond.attrib["type1"] = bond.attrib["class1"]
        bond.attrib["type2"] = bond.attrib["class2"]
        del bond.attrib["class1"]
        del bond.attrib["class2"]

        # Fix prefixes
        types = [bond.get("type1"), bond.get("type2")]
        if any(t in glycam_types for t in types):
            bond.set("type1", replacements[types[0]])
            bond.set("type2", replacements[types[1]])
        else:
            force.remove(bond)

    # Process <HarmonicAngleForce>

    force = root.find("HarmonicAngleForce")
    for angle in force.findall("Angle"):
        # Change attributes from class to type
        angle.attrib["type1"] = angle.attrib["class1"]
        angle.attrib["type2"] = angle.attrib["class2"]
        angle.attrib["type3"] = angle.attrib["class3"]
        del angle.attrib["class1"]
        del angle.attrib["class2"]
        del angle.attrib["class3"]

        # Fix prefixes
        types = [angle.get("type1"), angle.get("type2"), angle.get("type3")]
        if any(t in glycam_types for t in types):
            angle.set("type1", replacements[types[0]])
            angle.set("type2", replacements[types[1]])
            angle.set("type3", replacements[types[2]])
        else:
            force.remove(angle)

    # Process <PeriodicTorsionForce>

    force = root.find("PeriodicTorsionForce")
    for tag in ["Proper", "Improper"]:
        for torsion in force.findall(tag):
            # Change attributes from class to type and remove bad torsion
            torsion.attrib["type1"] = torsion.attrib["class1"]
            torsion.attrib["type2"] = torsion.attrib["class2"]
            torsion.attrib["type3"] = torsion.attrib["class3"]
            torsion.attrib["type4"] = torsion.attrib["class4"]
            del torsion.attrib["class1"]
            del torsion.attrib["class2"]
            del torsion.attrib["class3"]
            del torsion.attrib["class4"]
            if (
                torsion.attrib["type1"] == "NH"
                and torsion.attrib["type2"] == "Cg"
                and torsion.attrib["type3"] == "Cg"
                and torsion.attrib["type4"] == "Sm"
            ):
                force.remove(torsion)
                continue

            # Fix prefixes
            types = [
                torsion.get("type1"),
                torsion.get("type2"),
                torsion.get("type3"),
                torsion.get("type4"),
            ]
            if any(t in glycam_types for t in types):
                torsion.set("type1", replacements[types[0]])
                torsion.set("type2", replacements[types[1]])
                torsion.set("type3", replacements[types[2]])
                torsion.set("type4", replacements[types[3]])
            else:
                force.remove(torsion)

    # Process <NonbondedForce>

    force = root.find("NonbondedForce")
    for atom in force.findall("Atom"):
        # Change attributes from class to type
        atom.attrib["type"] = atom.attrib["class"]
        del atom.attrib["class"]

        # Fix prefixes
        type = atom.get("type")
        if type in glycam_types:
            atom.set("type", replacements[type])
        else:
            force.remove(atom)

    # Remove bad NH torsion and update the unscaled types in the script

    script = root.find("Script")
    text = script.text

    # Bad NH torsion
    text_split = text.split("\n")
    text_NH_removed = [line for line in text_split if "NH" not in line]
    text = "\n".join(text_NH_removed)

    # Unscaled types
    pattern = re.compile("unscaled_types = set\\((\\[(.*\n)*?.*?\\])\\)")
    match = pattern.search(text)
    types = eval(match.group(1))
    types = [(replacements[t[0]], replacements[t[1]], replacements[t[2]], replacements[t[3]]) for t in types]
    types = ",\n    ".join(str(t) for t in types)
    script.text = pattern.sub(f"unscaled_types = set([{types}])", text)

    # Add initialization script for setting up GlycamTemplateMatcher
    initialization_script = etree.SubElement(root, "InitializationScript")
    initialization_script.text = """
from openmm.app.internal import compiled

class GlycamTemplateMatcher(object):

  def __init__(self, glycam_residues):
    self.glycam_residues = glycam_residues

  def __call__(self, ff, residue, bondedToAtom, ignoreExternalBonds, ignoreExtraParticles):
    if residue.name in self.glycam_residues:
      template = ff._templates[residue.name]
      if compiled.matchResidueToTemplate(
        residue,
        template,
        bondedToAtom,
        ignoreExternalBonds,
        ignoreExtraParticles,
    ) is not None:
        return template

      # The residue doesn't actually match the template with the same name.  Try the terminal variants.

      if 'N'+residue.name in self.glycam_residues:
        template = ff._templates['N'+residue.name]
        if compiled.matchResidueToTemplate(
            residue,
            template,
            bondedToAtom,
            ignoreExternalBonds,
            ignoreExtraParticles,
        ) is not None:
          return template
      if 'C'+residue.name in self.glycam_residues:
        template = ff._templates['C'+residue.name]
        if compiled.matchResidueToTemplate(
            residue,
            template,
            bondedToAtom,
            ignoreExternalBonds,
            ignoreExtraParticles,
        ) is not None:
          return template
    return None

glycam_residues = set()
for residue in tree.getroot().find('Residues').findall('Residue'):
  glycam_residues.add(residue.get('name'))
self.registerTemplateMatcher(GlycamTemplateMatcher(glycam_residues))

"""

    tree.write(input_ffxml_path)


def validate_glyco_protein(
    ffxml_name,
    leaprc_name,
    supp_leaprc_name="oldff/leaprc.ff14SB",
    supp_ffxml_name="../openmmforcefields/ffxml/amber/protein.ff14SB.xml",
):
    modify_glycan_ffxml(ffxml_name)

    if verbose:
        print(f"Glycosylated protein energy validation for {ffxml_name}")
    top = "files/glycam/Glycoprotein_shortened.parm7"
    crd = "files/glycam/Glycoprotein_shortened.rst7"
    assert_energies_glyco_protein(top, crd, (supp_ffxml_name, ffxml_name))
    if verbose:
        print(f"Glycosylated protein energy validation for {ffxml_name} was successful!")


def validate_nucleic(
    ffxml_name,
    args,
):
    raise NotImplementedError()


class Logger:
    """
    Log energy discrepancies to a file.

    Parameters
    ----------
    log_filename : str
        Name of CSV file to write to

    """

    # logs testing energies into csv
    def __init__(self, log_filename=None):
        if log_filename:
            csvfile = open(log_filename, "w")
            fieldnames = [
                "ffxml_name",
                "data_type",
                "test_system",
                "units",
                "HarmonicBondForce",
                "HarmonicAngleForce",
                "PeriodicTorsionForce_dihedrals",
                "PeriodicTorsionForce_impropers",
                "NonbondedForce",
            ]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            self.csvfile = csvfile
            self.writer = writer
        else:
            self.csvfile = None
            self.writer = None

    def close(self):
        if self.csvfile:
            self.csvfile.close()

    def log(self, energies):
        if self.writer:
            self.writer.writerow(energies)


if __name__ == "__main__":
    main()
