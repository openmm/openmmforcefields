# CHARMM forcefield conversion tools for OpenMM

This directory contains files and scripts needed to convert the CHARMM forcefield to OpenMM `ffxml` files
Updated to July 2020 release from http://mackerell.umaryland.edu/charmm_ff.shtml#charmm

## Manifest
* `ffxml/` - converted OpenMM `ffxml` files for CHARMM36
* `tests/` - test systems for CHARMM36
* `charmm36.yaml` - yaml input file needed to drive bundling by `convert_charmm.py`
* `convert_charmm.py` - script to convert CHARMM `top` and `par` files to `ffxml`
* `test_charmm.py` - script to test CHARMM conversion after `convert_charmm.py` has been run
* `energy.py` - validate energies against CHARMM lite academic version inside a docker container (see below)

## Notes

Notes on files that were excluded from conversion:

* There are two glycolipid stream files with duplicate dihedrals with different values.
According to discussion with Alex MacKerell, the carb glycolipid file should be used so the lipid glycolipid stream file was excluded.
* `toppar_all36_prot_aldehydes.str` and `toppar_all36_na_modifications.str` have different values for the angle of atom types O CD CT2.
These files should not be used in the same system so both were excluded.  A new atom type is needed to correct this.
If needed, [CGenFF](https://cgenff.paramchem.org/) can be used for aldehydes or the user can convert these files at their own risk.

## Converting

Retrieve and unpack the CHARMM files
```
wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul20.tgz
tar zxf toppar_c36_jul20.tgz
```
Convert the solvent force fields
```
python convert_charmm.py --verbose --in files/waters.yaml
```
Convert non-solvent force fields
```
python convert_charmm.py --verbose --in files/charmm36.yaml
```

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
