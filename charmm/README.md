# CHARMM forcefield conversion tools for OpenMM

This directory contains files and scripts needed to convert the CHARMM forcefield to OpenMM `ffxml` files

## Manifest
* `toppar/` - CHARMM36 `toppar` files taken and unzipped from (http://mackerell.umaryland.edu/charmm_ff.html) Jan 2016.
* `ffxml/` - converted OpenMM `ffxml` files for CHARMM36
* `charmm36.yaml` - yaml input file needed to drive bundling by `convert_charmm.py`
* `convert_charmm.py` - script to convert CHARMM `top` and `par` files to `ffxml`

## Dependencies
* [`ParmEd`](https://github.com/parmed/parmed) - for parameter/topology conversion

## Notes

Notes on files that were excluded from conversion:

* There are two glycolipid stream files with duplicate dihedrals with different values.
According to discussion with Alex MacKerell, the carb glycolipid file should be used so the lipid glycolipid stream file was excluded.
* `toppar_all36_prot_aldehydes.str` and `toppar_all36_na_modifications.str` have different values for the angle of atom types O CD CT2.
These files should not be used in the same system so both were excluded.  A new atom type is needed to correct this.
If needed, [CGenFF](https://cgenff.paramchem.org/) can be used for aldehydes or the user can convert these files at their own risk.
