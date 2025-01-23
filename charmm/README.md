# CHARMM forcefield conversion tools for OpenMM

This directory contains files and scripts needed to convert the CHARMM forcefield to OpenMM `ffxml` files.

Updated to July 2024 release from <http://mackerell.umaryland.edu/charmm_ff.shtml#charmm>.

## Manifest
* `tests/` - test systems for CHARMM36
* `charmm36.yaml` - yaml input file needed to drive bundling by `convert_charmm.py`
* `convert_charmm.py` - script to convert CHARMM `top` and `par` files to `ffxml`
* `test_charmm.py` - script to test CHARMM conversion after `convert_charmm.py` has been run
* `energy.py` - validate energies against CHARMM lite academic version inside a docker container (see below)

## Notes

Notes on files or residues that were excluded from conversion due to
inconsistencies or other problems:

* `toppar/stream/na/toppar_all36_na_reactive_rna.str` is excluded due to a
  naming collision between atoms in residue DMPR (dimethylpropanamide) and patch
  DMPR (thio-substitution of dimyristoyl-D-glycero-1-phosphatidic acid).
* `toppar/stream/carb/toppar_all36_carb_glycolipid.str` is excluded due to
  parameter values for CG2R61-CG301 (bond from 6-membered aromatic ring carbon
  to aliphatic carbon with no hydrogens) that are inconsistent with the CGenFF
  values.
* `toppar/drude/drude_toppar_2023/toppar_drude_nucleic_acid_2017d.str` is
  excluded as ParmEd is unable to understand some of the entries in the file.
* Thulium(III) ion from `toppar/stream/misc/toppar_ions_won.str` (TM3P) is
  excluded because it collides with the identically-named 4'-methyl,3'-phosphate
  tetrahydrofuran residue.

If you need parameters from these files that were excluded, you may download the
CHARMM TOPPAR files, manually edit the offending entries, and regenerate the
OpenMM XML files at your own risk.

## Converting

Retrieve and unpack the CHARMM files
```
wget -O toppar.tgz http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
tar -xzf toppar.tgz
```
Convert the solvent force fields
```
python convert_charmm.py --verbose --in files/waters.yaml
```
Convert non-solvent force fields
```
python convert_charmm.py --verbose --in files/charmm36.yaml
```
Convert Drude force fields
```
(cd toppar/drude; tar -xzf drude_toppar_2023.tgz)
python convert_charmm.py --verbose --in files/drude2023.yaml
```

It is strongly advised to review the input conversion specification files (the
.yaml files) and all warning messages produced during conversion when updating
files to support a new version of the CHARMM force field parameters.
Idiosyncracies of each version of the force field parameters or a lack of
support from ParmEd may necessitate excluding certain parameter files.  If
warning messages indicating the presence of inconsistent parameter values are
generated, the TOPPAR files should be reviewed to locate the inconsistency.  You
may wish to exclude one or another set of parameters at this point to ensure
that the resulting force field is correct.

CHARMM docker image
===================

Prerequisites
-------------
* [Docker toolbox](https://www.docker.com/products/docker-toolbox)
* [CHARMM lite nonprofit/academic version](http://charmm.chemistry.harvard.edu/charmm_lite.php) (free) - downloaded as `charmm.tar.gz`

Building the docker image
-------------------------
After starting the Docker daemon, run
```
docker build -t omnia/charmm-lite:c40b1 .
```

Running CHARMM
--------------
The CHARMM executable is `/charmm/c40b1_gnu/exec/gnu/charmm`

To manually start the docker image (for testing purposes):
```
docker run -i -t omnia/charmm-lite:c40b1 /bin/bash

Running comparison
------------------
```
To run the comparison from this directory:
```
python energy.py dhfr
```
