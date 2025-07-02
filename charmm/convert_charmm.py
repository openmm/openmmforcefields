from parmed import openmm
from parmed.charmm import CharmmParameterSet
import argparse
import copy
import glob
import hashlib
import io
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
        help='path of the output directory. Default: "ffxml/" for yaml, "./" for leaprc',
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

        # compile residue and patch template names to exclude
        exclude_residue_names = source_files.get("exclude_residues", [])
        exclude_patch_names = source_files.get("exclude_patches", [])

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

        if verbose:
            print(f"Loading CHARMM parameter sets {charmm_files}...")
        params = CharmmParameterSet(*charmm_files)

        if exclude_residue_names:
            if verbose:
                print(f"Excluding residues: {exclude_residue_names}")
            for residue_name in exclude_residue_names:
                del params.residues[residue_name]
        if exclude_patch_names:
            if verbose:
                print(f"Excluding patches: {exclude_patch_names}")
            for patch_name in exclude_patch_names:
                del params.patches[patch_name]

        # Identical residue and patch names will cause problems, so ensure that
        # there are none present.
        name_collisions = set(params.residues) & set(params.patches)
        if name_collisions:
            raise ValueError(f"Collision between residue/patch names {sorted(name_collisions)}")

        # Ensure that no names contain hyphens, as this will confuse the
        # improper/anisotropy handling scripts.
        for templates in (params.residues.values(), params.patches.values()):
            for template in templates:
                if "-" in template.name:
                    raise ValueError(f"Forbidden character in template {template.name}")
                for atom in template.atoms:
                    if "-" in atom.name:
                        raise ValueError(f"Forbidden character in template {template.name} atom {atom.name}")

        # ParmEd does not handle CHARMM impropers properly, so save all of the
        # information about the impropers before discarding it from the
        # parameter set entirely.
        residue_impropers = {residue_name: residue._impr.copy() for residue_name, residue in params.residues.items()}
        patch_impropers = {patch_name: patch._impr.copy() for patch_name, patch in params.patches.items()}
        delete_impropers = {patch_name: patch.delete_impropers.copy() for patch_name, patch in params.patches.items()}
        periodic_improper_types = params.improper_periodic_types.copy()
        harmonic_improper_types = params.improper_types.copy()

        for residue in params.residues.values():
            residue._impr.clear()
        for patch in params.patches.values():
            patch._impr.clear()
        for patch in params.patches.values():
            patch.delete_impropers.clear()
        params.improper_periodic_types.clear()
        params.improper_types.clear()

        # ParmEd does not handle CHARMM anisotropies properly, so save all of
        # the information about them to deal with later.  The anisotropy
        # information will be removed from the Drude force entry in the force
        # field and applied in a script.
        residue_anisotropies = {
            residue_name: residue.anisotropies.copy() for residue_name, residue in params.residues.items()
        }
        patch_anisotropies = {patch_name: patch.anisotropies.copy() for patch_name, patch in params.patches.items()}
        delete_anisotropies = {
            patch_name: patch.delete_anisotropies.copy() for patch_name, patch in params.patches.items()
        }

        # For each entry, store the name of the XML file to contain the
        # residues.  Load CHARMM parameter sets to determine residues and
        # patches to split out.
        split_data = []
        if "split" in source_files:
            for split_spec in source_files["split"]:
                split_params = CharmmParameterSet(*split_spec["input"])
                split_fixes = split_spec.get("fixes", [])
                split_data.append(
                    (
                        os.path.join(ffxml_dir, split_spec["output"]),
                        set(split_params.residues),
                        set(split_params.patches),
                        split_fixes,
                    )
                )

        all_split_residue_names = set(
            residue_name for _, split_residue_names, _, _ in split_data for residue_name in split_residue_names
        )
        all_split_patch_names = set(
            patch_name for _, _, split_patch_names, _ in split_data for patch_name in split_patch_names
        )

        # Write a force field file out and read it back in.  This will contain
        # all residues and patches and will be larger than the "main" and
        # "split" force field files generated later.
        if verbose:
            print("Writing full force field file...")
        params_omm = openmm.OpenMMParameterSet.from_parameterset(params, unique_atom_types=True)
        params_omm.write(ffxml_filename, provenance=provenance, separate_ljforce=True)
        ffxml_tree = read_xml_file(ffxml_filename)

        # Prepare the main force field file, which will have some residues and
        # patches removed, then read it back in for verification.
        if verbose:
            print("Writing main force field file...")
        main_ffxml_tree = copy.deepcopy(ffxml_tree)
        main_residue_names, main_patch_names = strip_tree(
            main_ffxml_tree, all_split_residue_names, all_split_patch_names, False, set(), set()
        )
        write_improper_script(
            main_ffxml_tree,
            residue_impropers,
            patch_impropers,
            delete_impropers,
            periodic_improper_types,
            harmonic_improper_types,
        )
        write_anisotropy_script(
            main_ffxml_tree,
            residue_anisotropies,
            patch_anisotropies,
            delete_anisotropies,
        )

        apply_fixes(main_ffxml_tree, source_files.get("fixes", []))
        write_xml_file(main_ffxml_tree, ffxml_filename)
        app.ForceField(ffxml_filename)

        for split_ffxml_filename, split_residue_names, split_patch_names, split_fixes in split_data:
            # Prepare a split force field file, which will have some residues
            # and patches removed, and will contain only residues and patches,
            # then read it back in for verification.
            if verbose:
                print("Writing split force field file...")
            split_ffxml_tree = copy.deepcopy(ffxml_tree)
            strip_tree(
                split_ffxml_tree, split_residue_names, split_patch_names, True, main_residue_names, main_patch_names
            )
            apply_fixes(split_ffxml_tree, split_fixes)
            write_xml_file(split_ffxml_tree, split_ffxml_filename)
            app.ForceField(ffxml_filename, split_ffxml_filename)

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


# Allows edits to be made to the force field XML file in an automated way.  The
# target can use XPath notation to specify the location in the tree where the
# modification should take place.
def apply_fixes(ffxml_tree, fixes):
    for fix in fixes:
        action = fix["action"]
        target = fix.get("target", ".")

        if isinstance(target, list):
            target_elements = []
            for target_item in target:
                target_elements.extend(ffxml_tree.findall(target_item))
        else:
            target_elements = ffxml_tree.findall(target)

        if action == "append":
            content = fix["content"]
            if not isinstance(content, list):
                content = [content]

            for target_element in target_elements:
                for content_item in content:
                    target_element.append(build_xml_element(content_item))

        elif action == "remove":
            find_and_remove_xml_elements(ffxml_tree.getroot(), target_elements)

        else:
            raise ValueError(f"Unknown action {action!r}")


def find_and_remove_xml_elements(root, target_elements):
    root[:] = [child for child in root if child not in target_elements]
    for child in root:
        find_and_remove_xml_elements(child, target_elements)


def build_xml_element(data):
    element = etree.Element(data["tag"], data.get("attrib", {}))
    element.text = data.get("text", None)
    element.tail = data.get("tail", None)
    for child in data.get("children", []):
        element.append(build_xml_element(child))
    return element


def strip_tree(
    ffxml_tree, check_residue_names, check_patch_names, is_split, main_residue_names, main_patch_names=None
):
    # Helper function to decide whether or not to keep a template.  If a tree
    # for the main file is being processed, we are given all of the split names
    # and should keep names not in the set.  If a tree for the split file is
    # being processed, we are given the split names for that file and should
    # keep only the names in the set.
    def keep_template(template_name, check_template_names):
        return is_split == (template_name in check_template_names)

    # Helper function to extract a true patch name from a patch name in the XML
    # file (which will also contain a version number).
    def get_patch_name(patch_versioned_name):
        return patch_versioned_name.rsplit("_", maxsplit=1)[0]

    root = ffxml_tree.getroot()

    main_elements = []
    apply_patch_entries = []
    saved_residue_elements = {}
    saved_patch_elements = {}

    for main_element in root:
        if main_element.tag == "Info":
            # Do not touch <Info> tags.
            pass

        elif main_element.tag == "Residues":
            # Filter <Residues> tags.
            residue_elements = []

            for residue_element in main_element:
                residue_name = residue_element.attrib["name"]

                # We expect that ParmEd will have written information only in
                # <ApplyToResidue> tags, not <AllowPatch> tags.
                for item_element in residue_element:
                    if item_element.tag == "AllowPatch":
                        raise ValueError("expected only <ApplyToResidue>, not <AllowPatch>")

                if not keep_template(residue_name, check_residue_names):
                    continue

                saved_residue_elements[residue_name] = residue_element
                residue_elements.append(residue_element)

            if not residue_elements:
                continue
            main_element[:] = residue_elements

        elif main_element.tag == "Patches":
            # Filter <Patches> tags.
            patch_elements = []

            for patch_element in main_element:
                patch_versioned_name = patch_element.attrib["name"]

                # Save the information from all <ApplyToResidue> tags before
                # discarding them.
                item_elements = []
                for item_element in patch_element:
                    if item_element.tag == "ApplyToResidue":
                        apply_patch_entries.append((patch_versioned_name, item_element.attrib["name"]))
                        continue
                    item_elements.append(item_element)
                patch_element[:] = item_elements

                if not keep_template(get_patch_name(patch_versioned_name), check_patch_names):
                    continue

                saved_patch_elements[patch_versioned_name] = patch_element
                patch_elements.append(patch_element)

            if not patch_elements:
                continue
            main_element[:] = patch_elements

        elif main_element.tag == "DrudeForce":
            # Remove <DrudeForce> from split files; remove anisotropy
            # information from main files.
            if is_split:
                continue
            for particle_element in main_element:
                for anisotropy_attrib in ("type3", "type4", "type5", "aniso12", "aniso34"):
                    particle_element.attrib.pop(anisotropy_attrib, None)

        else:
            # Remove all other tags from split files.
            if is_split:
                continue

        main_elements.append(main_element)

    root[:] = main_elements

    # Add back <AllowPatch> and <ApplyToResidue> tags if possible.
    for patch_versioned_name, residue_name in apply_patch_entries:
        if residue_name in saved_residue_elements and (
            patch_versioned_name in saved_patch_elements or get_patch_name(patch_versioned_name) in main_patch_names
        ):
            etree.SubElement(saved_residue_elements[residue_name], "AllowPatch", name=patch_versioned_name)
        elif patch_versioned_name in saved_patch_elements and (
            residue_name in saved_residue_elements or residue_name in main_residue_names
        ):
            etree.SubElement(saved_patch_elements[patch_versioned_name], "ApplyToResidue", name=residue_name)

    return set(saved_residue_elements), set(
        get_patch_name(patch_versioned_name) for patch_versioned_name in saved_patch_elements
    )


def write_improper_script(
    ffxml_tree, residue_impropers, patch_impropers, delete_impropers, periodic_improper_types, harmonic_improper_types
):
    # Helper function to remove a prefix from an atom name and return whether
    # the atom is in the previous (-1) residue, the current (0) residue, or the
    # next (+1) residue, along with the stripped atom name.
    def preprocess_atom(atom):
        if atom.startswith("-"):
            return (-1, atom[1:])
        if atom.startswith("+"):
            return (1, atom[1:])
        return (0, atom)

    # Helper function to replace a wildcard "X" with None.
    def preprocess_atom_type(atom_type):
        return None if atom_type == "X" else atom_type

    theta_conversion = unit.degree.conversion_factor_to(unit.radian)
    k_conversion = unit.kilocalorie.conversion_factor_to(unit.kilojoule)

    # Prepares arguments for PeriodicTorsionForce.addTorsion() from a
    # DihedralType object.
    def prepare_periodic_parameters(improper_type):
        return (improper_type.per, improper_type.phase * theta_conversion, improper_type.phi_k * k_conversion)

    # Prepares arguments for CustomTorsionForce.addTorsion() from an
    # ImproperType object.
    def prepare_harmonic_parameters(improper_type):
        return ((improper_type.psi_k * k_conversion, improper_type.psi_eq * theta_conversion),)

    # Helper function to remove prefixes from atom names and remove entries for
    # templates with no impropers.
    def preprocess_template_impropers(template_impropers):
        return {
            residue: tuple(tuple(preprocess_atom(atom) for atom in improper) for improper in impropers)
            for residue, impropers in template_impropers.items()
            if impropers
        }

    # Helper function to prepare force field parameters and remove wildcard "X"
    # characters from improper type specifications.
    def preprocess_improper_types(improper_types, prepare_parameters):
        result = {}

        for atom_types, improper_type in improper_types.items():
            key = tuple(preprocess_atom_type(atom_type) for atom_type in atom_types)
            value = prepare_parameters(improper_type)

            # Add the forward and reverse versions of the improper for easier
            # lookup by the script.  Ensure that there are no duplicates.
            for result_key in (key, key[::-1]):
                if result_key in result:
                    if result[result_key] != value:
                        raise ValueError(f"conflict for improper {result_key}")
                else:
                    result[result_key] = value

        return result

    residue_impropers = preprocess_template_impropers(residue_impropers)
    patch_impropers = preprocess_template_impropers(patch_impropers)
    delete_impropers = preprocess_template_impropers(delete_impropers)

    residue_order = {residue: index for index, residue in enumerate(residue_impropers)}
    patch_order = {patch: index for index, patch in enumerate(patch_impropers)}

    periodic_improper_types = preprocess_improper_types(periodic_improper_types, prepare_periodic_parameters)
    harmonic_improper_types = preprocess_improper_types(harmonic_improper_types, prepare_harmonic_parameters)

    # Do not write a script to the file if there are no impropers present.
    if not (residue_impropers or patch_impropers):
        return

    with open("convert_charmm_improper_script.txt") as script_file:
        script_template = script_file.read()

    script = etree.SubElement(ffxml_tree.getroot(), "Script")
    script.text = script_template.format(
        residue_impropers=format_dict(residue_impropers),
        patch_impropers=format_dict(patch_impropers),
        delete_impropers=format_dict(delete_impropers),
        residue_order=format_dict(residue_order),
        patch_order=format_dict(patch_order),
        periodic_improper_types=format_dict(periodic_improper_types),
        harmonic_improper_types=format_dict(harmonic_improper_types),
    )


def write_anisotropy_script(ffxml_tree, residue_anisotropies, patch_anisotropies, delete_anisotropies):
    # Helper function to extract data from anisotropy objects.
    def preprocess_template_anisotropies(template_anisotropies):
        return {
            residue: {
                anisotropy.atom1.name: (
                    anisotropy.atom2.name,
                    anisotropy.atom3.name,
                    anisotropy.atom4.name,
                    anisotropy.a11,
                    anisotropy.a22,
                )
                for anisotropy in anisotropies
            }
            for residue, anisotropies in template_anisotropies.items()
            if anisotropies
        }

    residue_anisotropies = preprocess_template_anisotropies(residue_anisotropies)
    patch_anisotropies = preprocess_template_anisotropies(patch_anisotropies)
    delete_anisotropies = {
        residue: [anisotropy[0] for anisotropy in anisotropies]
        for residue, anisotropies in delete_anisotropies.items()
        if anisotropies
    }

    residue_order = {residue: index for index, residue in enumerate(residue_anisotropies)}
    patch_order = {patch: index for index, patch in enumerate(patch_anisotropies)}

    # Do not write a script to the file if there are no anisotropies present.
    if not (residue_anisotropies or patch_anisotropies):
        return

    with open("convert_charmm_anisotropy_script.txt") as script_file:
        script_template = script_file.read()

    script = etree.SubElement(ffxml_tree.getroot(), "Script")
    script.text = script_template.format(
        residue_anisotropies=format_dict(residue_anisotropies),
        patch_anisotropies=format_dict(patch_anisotropies),
        delete_anisotropies=format_dict(delete_anisotropies),
        residue_order=format_dict(residue_order),
        patch_order=format_dict(patch_order),
    )


def format_dict(d):
    return "\n".join(["{"] + [f"    {k!r}: {v!r}," for k, v in d.items()] + ["}"])


if __name__ == "__main__":
    main()
