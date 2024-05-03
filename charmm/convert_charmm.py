from parmed.charmm import CharmmParameterSet
from parmed import openmm
import glob
import yaml
import hashlib
import os
import openmm.app as app
import argparse
import csv


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

        # generate recommended combination for charmm36
        if verbose:
            print(f"Loading CHARMM parameter sets {charmm_files}...")
        params = CharmmParameterSet(*charmm_files)

        if len(exclude_residues) > 0:
            if verbose:
                print(f"Excluding residues: {exclude_residues}")
            for resname in exclude_residues:
                del params.residues[resname]

        if verbose:
            print("Converting parameters to OpenMM...")
        params_omm = openmm.OpenMMParameterSet.from_parameterset(params, unique_atom_types=True)

        # Set override level
        if "Override" in entry:
            override_level = int(entry["Override"])
            if verbose:
                print("Setting residues and patches to override level %d..." % override_level)
            for name, residue in params_omm.residues.items():
                residue.override_level = override_level
            for name, patch in params_omm.patches.items():
                patch.override_level = override_level

        if verbose:
            print("Writing parameter set and compatible patches. This may take several minutes...")
        params_omm.write(ffxml_filename, provenance=provenance, charmm_imp=True, separate_ljforce=True)

        # Try reading the forcefield back in to make sure it is valid
        if verbose:
            print("Verifying ffxml file integrity...")
        if "TestInclude" in entry:
            ffxml_include = entry["TestInclude"]
            forcefield = app.ForceField(ffxml_filename, *ffxml_include)
        else:
            forcefield = app.ForceField(ffxml_filename)

        if "Test" in entry:
            for filename in entry["Test"]:
                if verbose:
                    print(f"Testing with {filename} ...")
                pdbfile = app.PDBFile(filename)
                forcefield.createSystem(pdbfile.topology)

        if verbose:
            print("Verified.")


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


if __name__ == "__main__":
    main()
