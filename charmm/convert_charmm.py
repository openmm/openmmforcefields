from parmed import modeller, openmm
from parmed.charmm import CharmmParameterSet
import argparse
import glob
import hashlib
import io
import itertools
import math
import openmm.app as app
import openmm.unit as unit
import os
import xml.etree.ElementTree as etree
import yaml


def main():
    global verbose
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
    args = parser.parse_args()
    verbose = args.verbose

    convert_yaml(args.input, ffxml_dir=args.output_dir)


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

        # Loads the main CHARMM parameter set.  This is abstracted into a
        # function because we need a few copies that get mutated in different
        # ways, and ParmEd has issues with deepcopy(), so we just reload it.
        def load_params(return_impropers=False):
            if verbose:
                print(f"Loading CHARMM parameter sets {charmm_files}...")
            params = CharmmParameterSet(*charmm_files)

            if len(exclude_residues) > 0:
                if verbose:
                    print(f"Excluding residues: {exclude_residues}")
                for resname in exclude_residues:
                    del params.residues[resname]

            name_collisions = set(params.residues) & set(params.patches)
            if name_collisions:
                raise ValueError(f"Collision between residue/patch names {sorted(name_collision)}")

            # ParmEd does not handle CHARMM impropers properly, so save all of the
            # information about the impropers before discarding it from the
            # parameter set entirely
            residue_impropers = {
                residue_name: residue._impr.copy() for residue_name, residue in params.residues.items()
            }
            patch_impropers = {patch_name: patch._impr.copy() for patch_name, patch in params.patches.items()}
            periodic_improper_types = params.improper_periodic_types.copy()
            harmonic_improper_types = params.improper_types.copy()

            for residue in params.residues.values():
                residue._impr.clear()
            for patch in params.patches.values():
                patch._impr.clear()
            params.improper_periodic_types.clear()
            params.improper_types.clear()

            if return_impropers:
                return params, residue_impropers, patch_impropers, periodic_improper_types, harmonic_improper_types
            else:
                return params

        # Load CHARMM parameter sets to determine residues and patches to split
        # For each entry, store the name of the XML file to contain the residues
        # and patches split out, the residue names to split out, and the patch
        # names to split out
        split_data = []
        if "split" in source_files:
            for split_spec in source_files["split"]:
                split_params = CharmmParameterSet(*split_spec["input"])
                split_fixes = split_spec.get("fixes", [])
                split_data.append(
                    (
                        os.path.join(ffxml_dir, split_spec["output"]),
                        list(split_params.residues.keys()),
                        list(split_params.patches.keys()),
                        split_fixes,
                    )
                )

        # Convert everything together to OpenMM so we can figure out which
        # patches apply to which residues
        if verbose:
            print("Determining patch applicability...")
        params_omm = openmm.OpenMMParameterSet.from_parameterset(load_params(), unique_atom_types=True)
        valid_residues_for_patch, valid_patches_for_residue = params_omm._determine_valid_patch_combinations(
            params_omm._find_unused_residues()
        )

        # For the main file, remove all residues and patches that are supposed
        # to go into the split files, convert from a CHARMM to an OpenMM
        # parameter set, write it out, then read it back in to check for
        # validity.
        if verbose:
            print("Preparing to write main force field file...")

        params_main, residue_impropers, patch_impropers, periodic_improper_types, harmonic_improper_types = (
            load_params(return_impropers=True)
        )

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
        params_main_omm.write(ffxml_filename, provenance=provenance, separate_ljforce=True)

        if verbose:
            print("Applying fixes to main force field file...")
        ffxml_tree = read_xml_file(ffxml_filename)
        apply_fixes(ffxml_tree, source_files.get("fixes", []))
        apply_impropers(
            ffxml_tree,
            params_main_omm,
            residue_impropers,
            patch_impropers,
            periodic_improper_types,
            harmonic_improper_types,
            valid_residues_for_patch,
            valid_patches_for_residue,
        )
        write_xml_file(ffxml_tree, ffxml_filename)

        if verbose:
            print("Reading back main force field file...")
        main_forcefield = app.ForceField(ffxml_filename, *ffxml_include)

        # Include only those residues and patches that belong in the split files
        for split_ffxml_filename, split_residue_names, split_patch_names, split_fixes in split_data:
            if verbose:
                print("Preparing to write split force field file...")

            # To ensure that all appropriate parameters are present, write out
            # additional residues and patches applicable to the patches and
            # residues in the split file.  We will then remove the duplicated
            # entries from the written XML file.
            write_residue_names = set(split_residue_names) | {
                residue_name
                for patch_name in split_patch_names
                for residue_name in valid_residues_for_patch.get(patch_name, [])
            }
            write_patch_names = set(split_patch_names) | {
                patch_name
                for residue_name in split_residue_names
                for patch_name in valid_patches_for_residue.get(residue_name, [])
            }

            params_split = load_params()
            for residue_name in list(params_split.residues.keys()):
                if residue_name not in write_residue_names:
                    del params_split.residues[residue_name]
            for patch_name in list(params_split.patches.keys()):
                if patch_name not in write_patch_names:
                    del params_split.patches[patch_name]

            params_split_omm = openmm.OpenMMParameterSet.from_parameterset(params_split, unique_atom_types=True)

            # These should override residues and patches in the main file
            for residue in params_split_omm.residues.values():
                residue.override_level = override_level + 1
            for patch in params_split_omm.patches.values():
                patch.override_level = override_level + 1

            # Use charmm_imp=False: see prepare_write() below.
            if verbose:
                print("Writing split force field file...")
            params_split_omm.write(split_ffxml_filename, provenance=provenance, separate_ljforce=True)

            if verbose:
                print("Applying fixes to split force field file...")
            split_ffxml_tree = read_xml_file(split_ffxml_filename)
            apply_fixes(split_ffxml_tree, split_fixes)
            apply_impropers(
                split_ffxml_tree,
                params_split_omm,
                residue_impropers,
                patch_impropers,
                periodic_improper_types,
                harmonic_improper_types,
                valid_residues_for_patch,
                valid_patches_for_residue,
            )
            clean_split(ffxml_tree, split_ffxml_tree)
            write_xml_file(split_ffxml_tree, split_ffxml_filename)

            if verbose:
                print("Reading back split force field file...")
            app.ForceField(ffxml_filename, split_ffxml_filename, *ffxml_include)

        if "Test" in entry:
            for filename in entry["Test"]:
                if verbose:
                    print(f"Testing with {filename} ...")
                pdbfile = app.PDBFile(filename)
                main_forcefield.createSystem(pdbfile.topology)

        if verbose:
            print("Done.")


def read_xml_file(xml_filename):
    return etree.parse(xml_filename)


def write_xml_file(tree, xml_filename):
    tree = etree.parse(
        io.StringIO(etree.canonicalize(etree.tostring(tree.getroot(), encoding="unicode"), strip_text=True))
    )
    etree.indent(tree)
    tree.write(xml_filename, encoding="unicode")


def apply_fixes(tree, fixes):
    # Allows edits to be made to the force field XML file in an automated way.
    # The target can use XPath notation to specify the location in the tree
    # where the modification should take place.
    for fix in fixes:
        action = fix["action"]
        for target_element in tree.findall(fix.get("target", ".")):
            if action == "append":
                target_element.append(build_xml_element(fix["content"]))
            else:
                raise ValueError(f"Unknown action {action!r}")


def build_xml_element(data):
    element = etree.Element(data["tag"], data.get("attrib", {}))
    element.text = data.get("text", None)
    element.tail = data.get("tail", None)
    for child in data.get("children", []):
        element.append(build_xml_element(child))
    return element


def apply_impropers(
    tree,
    params_omm,
    residue_impropers,
    patch_impropers,
    periodic_improper_types,
    harmonic_improper_types,
    valid_residues_for_patch,
    valid_patches_for_residue,
):
    k_conversion = unit.kilocalorie.conversion_factor_to(unit.kilojoule)
    theta0_conversion = unit.degree.conversion_factor_to(unit.radian)

    periodic_improper_element = etree.Element("PeriodicTorsionForce", ordering="charmm")
    harmonic_improper_element = etree.Element(
        "CustomTorsionForce",
        ordering="charmm",
        energy=f"k*min(dtheta, 2*pi - dtheta)^2; dtheta = abs(theta - theta0); pi = {math.pi}",
    )

    # Find the residues that ParmEd actually wrote
    residue_names = []
    for residue_element in tree.findall("./Residues/Residue"):
        residue_names.append(residue_element.attrib["name"])
    residue_name_set = set(residue_names)

    # Find the patches that ParmEd actually wrote
    patch_names = []
    for patch_element in tree.findall("./Patches/Patch"):
        patch_names.append(patch_element.attrib["name"])
    patch_name_set = set(patch_names)

    # Convert valid_residues_for_patch, a dictionary from patch names to lists
    # of residue names, into patch_compatible_residues, a dictionary from patch
    # names to residue objects.  Filter by only the residues and patches that
    # ParmEd actually wrote to the XML file.
    patch_compatible_residues = {
        patch_name: [
            params_omm.residues[residue_name] for residue_name in residue_names if residue_name in residue_name_set
        ]
        for patch_name, residue_names in valid_residues_for_patch.items()
        if patch_name in patch_name_set
    }

    # Some patches might be applied to residues where atoms in a patch improper
    # normally found in the patched residue actually end up belonging to a
    # different patch.  Find the patches compatible with the residues compatible
    # with each patch.
    patch_compatible_patches = {
        patch_name: [
            params_omm.patches[patch_name]
            for patch_name in patch_names
            if patch_name
            in set(
                compatible_patch_name
                for compatible_residue in compatible_residues
                for compatible_patch_name in valid_patches_for_residue.get(compatible_residue.name, [])
            )
        ]
        for patch_name, compatible_residues in patch_compatible_residues.items()
    }

    # Look up all of the atom names in adjacent residues referred to by
    # impropers
    adjacent_names = set()
    for residue_name in residue_names:
        for atom_names in residue_impropers[residue_name]:
            for atom_name in atom_names:
                if is_adjacent_name(atom_name):
                    adjacent_names.add(get_adjacent_name(atom_name))
    for patch_name in patch_names:
        for atom_names in patch_impropers[patch_name]:
            for atom_name in atom_names:
                if is_adjacent_name(atom_name):
                    adjacent_names.add(get_adjacent_name(atom_name))

    # Look up all of the atom types that could correspond to these names
    adjacent_types = dict()
    for adjacent_name in adjacent_names:
        adjacent_types_for_name = set()
        for residue_name in residue_names:
            residue = params_omm.residues[residue_name]
            if adjacent_name in residue.map:
                adjacent_types_for_name.add(residue.map[adjacent_name].type)
        for patch_name in patch_names:
            patch = params_omm.patches[patch_name]
            if adjacent_name in patch.map:
                adjacent_types_for_name.add(patch.map[adjacent_name].type)
        adjacent_types[adjacent_name] = adjacent_types_for_name

    def get_periodic_attributes(improper_data):
        return {
            "periodicity1": str(improper_data.per),
            "k1": str(improper_data.phi_k * k_conversion),
            "phase1": str(improper_data.phase * theta0_conversion),
        }

    def get_harmonic_attributes(improper_data):
        return {"k": str(improper_data.psi_k * k_conversion), "theta0": str(improper_data.psi_eq * theta0_conversion)}

    def process_impropers(residue_or_patch, impropers):
        for atom_names in impropers:
            # Handle periodic impropers
            for type_names, improper_data in periodic_improper_types.items():
                force_attributes = get_periodic_attributes(improper_data)
                if improper_matches(
                    residue_or_patch,
                    atom_names,
                    type_names,
                    adjacent_types,
                    patch_compatible_residues,
                    patch_compatible_patches,
                ):
                    for match_attributes in expand_names(
                        residue_or_patch,
                        atom_names,
                        adjacent_types,
                        patch_compatible_residues,
                        patch_compatible_patches,
                    ):
                        etree.SubElement(periodic_improper_element, "Improper", **match_attributes, **force_attributes)

            # Handle harmonic impropers
            for type_names, improper_data in harmonic_improper_types.items():
                force_attributes = get_harmonic_attributes(improper_data)
                if improper_matches(
                    residue_or_patch,
                    atom_names,
                    type_names,
                    adjacent_types,
                    patch_compatible_residues,
                    patch_compatible_patches,
                ):
                    for match_attributes in expand_names(
                        residue_or_patch,
                        atom_names,
                        adjacent_types,
                        patch_compatible_residues,
                        patch_compatible_patches,
                    ):
                        etree.SubElement(harmonic_improper_element, "Improper", **match_attributes, **force_attributes)

    # Write residue and patch impropers
    for residue_name in residue_names:
        process_impropers(params_omm.residues[residue_name], residue_impropers[residue_name])
    for patch_name in patch_names:
        process_impropers(params_omm.patches[patch_name], patch_impropers[patch_name])

    # Add force elements to the tree only if they have children
    root = tree.getroot()
    if len(periodic_improper_element):
        root.append(periodic_improper_element)
    if len(harmonic_improper_element):
        harmonic_improper_element[:0] = [
            etree.Element("PerTorsionParameter", name="k"),
            etree.Element("PerTorsionParameter", name="theta0"),
        ]
        root.append(harmonic_improper_element)


def is_adjacent_name(atom_name):
    return atom_name.startswith("-") or atom_name.startswith("+")


def get_adjacent_name(atom_name):
    return atom_name[1:]


def improper_matches(
    residue_or_patch, atom_names, type_names, adjacent_types, patch_compatible_residues, patch_compatible_patches
):
    def name_type_match(atom_name, type_name):
        if is_adjacent_name(atom_name):
            return type_name in adjacent_types[get_adjacent_name(atom_name)]
        else:
            if atom_name in residue_or_patch.map:
                return type_name == residue_or_patch.map[atom_name].type
            elif isinstance(residue_or_patch, modeller.PatchTemplate):
                return any(
                    atom_name in compatible_residue_or_patch.map
                    and type_name == compatible_residue_or_patch.map[atom_name].type
                    for patch_compatible in (patch_compatible_residues, patch_compatible_patches)
                    for compatible_residue_or_patch in patch_compatible[residue_or_patch.name]
                )
            else:
                raise KeyError(atom_name)

    return any(
        all(
            type_name == "X" or name_type_match(atom_name, type_name)
            for atom_name, type_name in zip(atom_names[::order], type_names)
        )
        for order in (1, -1)
    )


def expand_names(residue_or_patch, atom_names, adjacent_types, patch_compatible_residues, patch_compatible_patches):
    def expand_name(atom_name, compatible_residue_or_patch):
        if is_adjacent_name(atom_name):
            for type_name in adjacent_types[get_adjacent_name(atom_name)]:
                yield "class", type_name
        else:
            if atom_name in residue_or_patch.map:
                # Atom in this residue or patch
                yield "type", f"{residue_or_patch.name}-{atom_name}"
            elif compatible_residue_or_patch is not None and atom_name in compatible_residue_or_patch.map:
                # Atom in the compatible residue or patch of this patch
                yield "type", f"{compatible_residue_or_patch.name}-{atom_name}"

    if isinstance(residue_or_patch, modeller.PatchTemplate):
        compatible_residues_and_patches = (
            patch_compatible_residues[residue_or_patch.name] + patch_compatible_patches[residue_or_patch.name]
        )
    else:
        compatible_residues_and_patches = [None]

    # If residue_or_patch is a residue, compatible_residues_and_patches is [None] and we go
    # generate duplicates, which are screened out.
    # patch contains impropers not involving atoms in the residue or patch, this will
    # patch, we go through for each residue or patch compatible with the patch.  If the
    # through this outer loop once, considering only the residue.  If it is a
    yielded = set()
    for compatible_residue_or_patch in compatible_residues_and_patches:
        for specification in itertools.product(
            *(expand_name(atom_name, compatible_residue_or_patch) for atom_name in atom_names)
        ):
            if specification in yielded:
                continue
            yielded.add(specification)
            yield {f"{key_kind}{index + 1}": value for index, (key_kind, value) in enumerate(specification)}


def element_to_key(element):
    return etree.canonicalize(etree.tostring(element, encoding="unicode"), strip_text=True)


def clean_split(main_ffxml_tree, split_ffxml_tree):
    # Find the names of the residues and patches in the main file.
    main_residues = set()
    for residue_element in main_ffxml_tree.findall("./Residues/Residue"):
        main_residues.add(residue_element.attrib["name"])
    main_patches = set()
    for patch_element in main_ffxml_tree.findall("./Patches/Patch"):
        main_patches.add(patch_element.attrib["name"])

    # Delete residues that are already in the main file.  If a residue is not in
    # the main file, leave it in the split file; if it is in the main file, save
    # its <AllowPatch> tags corresponding to patches in the split file before
    # deleting it.
    patches_apply_to_residue = {}
    for residues_element in split_ffxml_tree.findall("./Residues"):
        new_residue_elements = []
        for residue_element in residues_element:
            residue_name = residue_element.attrib["name"]
            if residue_name in main_residues:
                for allow_patch_element in residue_element.findall("./AllowPatch"):
                    patch_name = allow_patch_element.attrib["name"]
                    if patch_name not in main_patches:
                        patches_apply_to_residue.setdefault(patch_name, []).append(residue_name)
            else:
                new_residue_elements.append(residue_element)
        residues_element[:] = new_residue_elements

    # Delete patches that are already in the main file.  For patches retained in
    # the split file, add <ApplyToResidue> tags for any relevant residues in the
    # main file.
    for patches_element in split_ffxml_tree.findall("./Patches"):
        new_patch_elements = []
        for patch_element in patches_element:
            patch_name = patch_element.attrib["name"]
            if patch_name not in main_patches:
                for applicable_residue in patches_apply_to_residue.get(patch_name, []):
                    etree.SubElement(patch_element, "ApplyToResidue", name=applicable_residue)
                new_patch_elements.append(patch_element)
        patches_element[:] = new_patch_elements

    def improper_element_to_key(improper_element):
        return tuple(sorted(improper_element.attrib.items()))

    # Find impropers already in the main file.  This assumes that there is only
    # one CustomTorsionForce that the impropers are a part of.  Impropers in a
    # PeriodicTorsionForce are assumed to be handled by OpenMM.
    main_impropers = set()
    for improper_element in main_ffxml_tree.findall("./CustomTorsionForce/Improper"):
        main_impropers.add(improper_element_to_key(improper_element))

    # Delete associated impropers in the split file.
    for impropers_element in split_ffxml_tree.findall("./CustomTorsionForce"):
        impropers_element[:] = [
            element
            for element in impropers_element
            if improper_element.tag != "Improper" or improper_element_to_key(improper_element) in main_impropers
        ]


if __name__ == "__main__":
    main()
