import subprocess

MOLECULE_COUNT = 855
SYSTEM_NAME = "WATBOX"

SPECIFICATIONS = [
    ("waterbox-3-site-tip3p.psf", "toppar_water_ions.str", "TIP3"),
    ("waterbox-3-site-tip3p-pme-b.psf", "non_charmm/toppar_water_ions_tip3p_pme_b.str", "TP3B"),
    ("waterbox-3-site-tip3p-pme-f.psf", "non_charmm/toppar_water_ions_tip3p_pme_f.str", "TP3F"),
    ("waterbox-3-site-spc.psf", "non_charmm/toppar_water_ions_spc.str", "SPC"),
    ("waterbox-3-site-spc-e.psf", "non_charmm/toppar_water_ions_spc_e.str", "SPCE"),
    ("waterbox-4-site-tip4p.psf", "non_charmm/toppar_water_ions_tip4p.str", "TIP4"),
    ("waterbox-4-site-tip4p-ew.psf", "non_charmm/toppar_water_ions_tip4p_ew.str", "TP4E"),
    ("waterbox-4-site-tip4p-2005.psf", "non_charmm/toppar_water_ions_tip4p_2005.str", "TP45"),
    ("waterbox-5-site-tip5p.psf", "non_charmm/toppar_water_ions_tip5p.str", "TIP5"),
    ("waterbox-5-site-tip5p-ew.psf", "non_charmm/toppar_water_ions_tip5p_ew.str", "TP5E"),
]


def main():
    for psf_path, str_path, residue_name in SPECIFICATIONS:
        input_lines = [
            "* Test",
            "*",
            f"stream ../../toppar/{str_path}",
            f"read sequence {residue_name} {MOLECULE_COUNT}",
            f"generate {SYSTEM_NAME} noangle nodihedral",
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
