from parmed import openmm
from parmed.charmm import CharmmParameterSet
import argparse
import copy
import csv
import glob
import hashlib
import io
import openmm.app as app
import os
import xml.etree.ElementTree as etree
import yaml


def main():
    global verbose
    global no_log
    global logger
    # Set up parser
    parser = argparse.ArgumentParser(description="CHARMM --> OpenMM forcefield conversion script")
    parser.add_argument(
        "--input", "-i", default="files/waters.yaml", help='path of the input file. Default: "files/waters.yaml"'
    )
    parser.add_argument(
        "--output-dir",
        "-od",
        help="path of the output directory. " 'Default: "ffxml/" for yaml, "./" for leaprc',
        default="ffxml/",
    )
    parser.add_argument("--verbose", "-v", action="store_true", help="turns verbosity on")
    parser.add_argument("--no-log", action="store_true", help="turns logging of energies to log.csv off")
    args = parser.parse_args()
    verbose = args.verbose
    no_log = args.no_log

    if not no_log:
        logger = Logger("log.csv")
    convert_yaml(args.input, ffxml_dir=args.output_dir)
    if not no_log:
        logger.close()


def convert_yaml(yaml_filename, ffxml_dir):
    # Read YAML
    data = yaml.safe_load(open(yaml_filename))
    source_pack = data[0]["sourcePackage"]
    source_pack_ver = data[0]["sourcePackageVersion"]

    for entry in data[1:]:
        ffxml_filename = os.path.join(ffxml_dir, entry["Destination"])
        print(f"Generating {ffxml_filename}")
        charmm_references = entry["References"]
        source_files = entry["Source"]

        # files that should be excluded from conversion.
        exclude_files = set()
        if ("exclude" in source_files) and (source_files["exclude"] is not None):
            exclude_files = set(source_files["exclude"])

        # charmm36 main top and par files
        charmm_files = list()
        if ("include" in source_files) and (source_files["include"] is not None):
            charmm_files = source_files["include"]

        # add stream files
        if ("stream" in source_files) and (source_files["stream"] is not None):
            for files in source_files["stream"]:
                charmm_files.extend(glob.glob(files))

        # compile residue template names to exclude
        exclude_residues = list()
        if "exclude_residues" in source_files:
            for resname in source_files["exclude_residues"]:
                exclude_residues.append(resname)

        # exclude files from conversion, maintaining deterministic order
        for filename in exclude_files:
            try:
                charmm_files.remove(filename)
            except Exception:
                print(f'Specified excluded file "{filename}" does not appear in list of files')

        provenance = dict()
        source = provenance["Source"] = []
        for fi in charmm_files:
            source.append(dict())
            source[-1]["Source"] = fi
            md5 = hashlib.md5()
            with open(fi, "rb") as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]["md5hash"] = md5
            source[-1]["sourcePackage"] = source_pack
            source[-1]["sourcePackageVersion"] = source_pack_ver

        references = provenance["Reference"] = []
        for ff in charmm_references:
            for cite in charmm_references[ff]:
                references.append(dict())
                if isinstance(cite, dict):
                    for key in cite.keys():
                        citation = cite[key]
                        references[-1]["Reference"] = citation
                        references[-1]["forcefield"] = ff
                        references[-1]["type"] = key
                else:
                    citation = cite
                    references[-1]["Reference"] = citation
                    references[-1]["forcefield"] = ff
 
        if "Override" in entry:
            override_level = int(entry["Override"])
            if verbose:
                print("Using override level %d..." % override_level)
        else:
            override_level = 0

        if "TestInclude" in entry:
            ffxml_include = entry["TestInclude"]
        else:
            ffxml_include = []
       
        # Load the main CHARMM parameter set
        if verbose:
            print(f"Loading CHARMM parameter sets {charmm_files}...")
        params = CharmmParameterSet(*charmm_files)
        
        if len(exclude_residues) > 0:
            if verbose:
                print(f"Excluding residues: {exclude_residues}")
            for resname in exclude_residues:
                del params.residues[resname]

        # Load CHARMM parameter sets to determine residues and patches to split
        # For each entry, store the name of the XML file containing the residues
        # and patches split out, the residue names to split out, and the patch
        # names to split out
        split_data = []
        if "split" in source_files:
            for split_spec in source_files["split"]:
                split_params = CharmmParameterSet(*split_spec["input"])
                split_fixes = split_spec.get("fixes", [])
                split_data.append((os.path.join(ffxml_dir, split_spec["output"]), list(split_params.residues.keys()), list(split_params.patches.keys()), split_fixes))

        # Convert everything together to OpenMM so we can figure out which
        # patches apply to which residues
        if verbose:
            print("Determining patch applicability...")
        params_omm = openmm.OpenMMParameterSet.from_parameterset(params, unique_atom_types=True)
        valid_residues_for_patch, valid_patches_for_residue = params_omm._determine_valid_patch_combinations(params_omm._find_unused_residues())

        # For the main file, remove all residues and patches that are supposed
        # to go into the split files, convert from a CHARMM to an OpenMM
        # parameter set, write it out, then read it back in to check for
        # validity.
        if verbose:
            print("Preparing to write main force field file...")

        params_main = copy.deepcopy(params)

        for split_ffxml_filename, split_residue_names, split_patch_names, split_fixes in split_data:
            for split_residue_name in split_residue_names:
                if split_residue_name in params_main.residues:
                    del params_main.residues[split_residue_name]
            for split_patch_name in split_patch_names:
                if split_patch_name in params_main.patches:
                    del params_main.patches[split_patch_name]

        params_main_omm = openmm.OpenMMParameterSet.from_parameterset(params_main, unique_atom_types=True)

        for residue in params_main_omm.residues.values():
            residue.override_level = override_level
        for patch in params_main_omm.patches.values():
            patch.override_level = override_level

        if verbose:
            print("Writing main force field file...")
        params_main_omm.write(ffxml_filename, provenance=provenance, charmm_imp=True, separate_ljforce=True)
        apply_fixes(ffxml_filename, source_files.get("fixes", []))

        if verbose:
            print("Reading back main force field file...")
        main_forcefield = app.ForceField(ffxml_filename, *ffxml_include)

        # For the split files, include only those residues and patches that
        # belong, plus any residues from the main file that have patches in the
        # split file applying to them.  ParmEd doesn't give us access to control
        # writing <ApplyToResidue>, unfortunately.
        for split_ffxml_filename, split_residue_names, split_patch_names, split_fixes in split_data:
            if verbose:
                print("Preparing to write split force field file...")

            write_residue_names = set(split_residue_names) | {residue_name for patch_name in split_patch_names for residue_name in valid_residues_for_patch.get(patch_name, ())}
            write_patch_names = set(split_patch_names)

            params_split = copy.deepcopy(params)
            for residue_name in list(params_split.residues.keys()):
                if residue_name not in write_residue_names:
                    del params_split.residues[residue_name]
            for patch_name in list(params_split.patches.keys()):
                if patch_name not in write_patch_names:
                    del params_split.patches[patch_name]

            params_split_omm = openmm.OpenMMParameterSet.from_parameterset(params_split, unique_atom_types=True)

            # These should override the residues and patches in the main file
            for residue in params_split_omm.residues.values():
                residue.override_level = override_level + 1
            for patch in params_split_omm.patches.values():
                patch.override_level = override_level + 1

            if verbose:
                print("Writing split force field file...")
            params_split_omm.write(split_ffxml_filename, provenance=provenance, charmm_imp=True, separate_ljforce=True)
            apply_fixes(split_ffxml_filename, split_fixes)

            if verbose:
                print("Reading back split force field file...")
            split_forcefield = app.ForceField(ffxml_filename, split_ffxml_filename, *ffxml_include)

        if "Test" in entry:
            for filename in entry["Test"]:
                if verbose:
                    print(f"Testing with {filename} ...")
                pdbfile = app.PDBFile(filename)
                main_forcefield.createSystem(pdbfile.topology)

        if verbose:
            print("Done.")

class Logger:
    # logs testing energies into csv
    def __init__(self, log_file):
        csvfile = open(log_file, "w")
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

    def close(self):
        self.csvfile.close()

    def log(self, energies):
        self.writer.writerow(energies)

def apply_fixes(xml_filename, fixes):
    # Allows things to be done to the force field XML file in an automated way.
    # Right now the only action is "append" to add a new element.  The target
    # can use XPath notation to specify the location in the tree where the
    # modification should take place.
    tree = etree.parse(xml_filename)
    for fix in fixes:
        action = fix["action"]
        for target_element in tree.findall(fix.get("target", ".")):
            if action == "append":
                target_element.append(build_xml_element(fix["content"]))
            else:
                raise ValueError(f"Unknown action {action!r}")

    # Pretty print the XML file
    tree = etree.parse(io.StringIO(etree.canonicalize(etree.tostring(tree.getroot(), encoding="unicode"), strip_text=True)))
    etree.indent(tree)
    tree.write(xml_filename, encoding="unicode")

def build_xml_element(data):
    element = etree.Element(data["tag"], data.get("attrib", {}))
    element.text = data.get("text", None)
    element.tail = data.get("tail", None)
    for child in data.get("children", []):
        element.append(build_xml_element(child))
    return element

if __name__ == "__main__":
    main()
