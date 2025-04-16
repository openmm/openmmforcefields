import subprocess

MOLECULE_COUNT = 855
SYSTEM_NAME = "WATBOX"

SPECIFICATIONS = [
    ("waterbox-3-site-tip3p.psf", "toppar_water_ions.str", "TIP3", False),
    ("waterbox-3-site-tip3p-pme-b.psf", "toppar_water_ions_tip3p_pme_b.str", "TP3B", False),
    ("waterbox-3-site-tip3p-pme-f.psf", "toppar_water_ions_tip3p_pme_f.str", "TP3F", False),
    ("waterbox-3-site-spc.psf", "toppar_water_ions_spc.str", "SPC", False),
    ("waterbox-3-site-spc-e.psf", "toppar_water_ions_spc_e.str", "SPCE", False),
    ("waterbox-4-site-tip4p.psf", "toppar_water_ions_tip4p.str", "TIP4", False),
    ("waterbox-4-site-tip4p-ew.psf", "toppar_water_ions_tip4p_ew.str", "TP4E", False),
    ("waterbox-4-site-tip4p-2005.psf", "toppar_water_ions_tip4p_2005.str", "TP45", False),
    ("waterbox-5-site-tip5p.psf", "toppar_water_ions_tip5p.str", "TIP5", False),
    ("waterbox-5-site-tip5p-ew.psf", "toppar_water_ions_tip5p_ew.str", "TP5E", False),
    ("waterbox-drude.psf", "toppar_drude_main_protein_2023a.str", "SWM4", True),
]


def main():
    for psf_path, str_path, residue_name, is_drude in SPECIFICATIONS:
        generate_extra = " drude" if is_drude else ""
        input_lines = [
            "* Test",
            "*",
            f"stream ../toppar/{str_path}",
            f"read sequence {residue_name} {MOLECULE_COUNT}",
            f"generate {SYSTEM_NAME} noangle nodihedral{generate_extra}",
            f"write psf card name {psf_path}",
            "* Water box",
            f"* Input file:  {str_path}",
            f"* Output file: {psf_path}",
            "*",
            "stop",
        ]

        subprocess.run(["charmm"], input="\n".join(input_lines).encode(), check=True)


if __name__ == "__main__":
    main()
