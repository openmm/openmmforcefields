import parmed

type_map = {
    "O": "OT",
    "OH2": "OT",
    "H1": "HT",
    "H2": "HT",
    "LP": "LP",
    "OM": "LP",
    "LP1": "LP",
    "LP2": "LP",
}


def rename_atoms(structure):
    for atom in structure.atoms:
        atom.type = type_map[atom.name]


def delete_lp_bonds(structure):
    indices_to_delete = []
    for index, bond in enumerate(structure.bonds):
        if (bond.atom1.type == "LP" and bond.atom2.type == "HT") or (
            bond.atom2.type == "LP" and bond.atom1.type == "HT"
        ):
            indices_to_delete.append(index)
    for index in indices_to_delete[::-1]:
        structure.bonds[index].delete()
        del structure.bonds[index]


def assign_charges(structure, parameter_set, residue_name):
    charges = {}
    for atom in parameter_set.residues[residue_name].atoms:
        if atom.type in charges:
            if charges[atom.type] != atom.charge:
                raise ValueError("inconsistent charges")
        else:
            charges[atom.type] = atom.charge
    for atom in structure.atoms:
        atom.charge = charges[atom.type]


def add_angles(structure):
    for residue in structure.residues:
        if not (residue.atoms[0].type == "OT" and residue.atoms[1].type == residue.atoms[2].type == "HT"):
            raise ValueError("unknown residue")
        structure.angles.append(parmed.Angle(residue.atoms[1], residue.atoms[0], residue.atoms[2]))


for site_count, str_path, residue_name, psf_path in [
    (3, "toppar_water_ions.str", "TIP3", "waterbox-3-site-tip3p.psf"),
    (3, "non_charmm/toppar_water_ions_spc.str", "SPC", "waterbox-3-site-spc.psf"),
    (3, "non_charmm/toppar_water_ions_spc_e.str", "SPCE", "waterbox-3-site-spc-e.psf"),
    (3, "non_charmm/toppar_water_ions_tip3p_pme_b.str", "TP3B", "waterbox-3-site-tip3p-pme-b.psf"),
    (3, "non_charmm/toppar_water_ions_tip3p_pme_f.str", "TP3F", "waterbox-3-site-tip3p-pme-f.psf"),
    # TODO: 4- and 5-site models
]:
    structure = parmed.load_file(f"waterbox-{site_count}-site.pdb")
    rename_atoms(structure)
    delete_lp_bonds(structure)

    psf = parmed.charmm.CharmmPsfFile.from_structure(structure)
    parameter_set = parmed.charmm.CharmmParameterSet(f"../../toppar/{str_path}")

    psf.load_parameters(parameter_set)
    assign_charges(structure, parameter_set, residue_name)
    add_angles(structure)

    parmed.formats.PSFFile.write(psf, psf_path)
