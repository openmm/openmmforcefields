# CHARMM forcefield conversion tools for OpenMM

This directory contains files and scripts needed to convert the CHARMM forcefield to OpenMM `ffxml` files.

Updated to July 2024 release from <http://mackerell.umaryland.edu/charmm_ff.shtml#charmm>.

## Manifest

* `convert_charmm.py`: Converts from CHARMM to FFXML force field format.
* `files`: Specification files for the conversions.
* `test_charmm.py`: Tests energies and forces from OpenMM FFXML vs. OpenMM
  reading CHARMM files vs. CHARMM itself.
* `tests`: Input files for test cases.

## Notes

Some files and residues are excluded from conversion due to inconsistencies or
other problems:

* `toppar/stream/carb/toppar_all36_carb_glycolipid.str` is excluded due to
  parameter values for CG2R61-CG301 (bond from 6-membered aromatic ring carbon
  to aliphatic carbon with no hydrogens) that are inconsistent with the CGenFF
  values.
* Residue DMPR (dimethylpropanamide) is excluded due to a collision with patch
  DMPR (thio-substitution of DMPA).
* Patch CNH3 (CYTP to protonated cytosine) is excluded because it is incorrectly
  marked as applicable to other residues, causing problems.
* `toppar/drude/drude_toppar_2023/toppar_drude_nucleic_acid_2017d.str` is
  excluded as ParmEd is unable to understand some of the entries in the file.
* Drude patches CTES and NTES are excluded, preventing parameterization of
  single amino acids, because inclusion of these patches causes patching
  ambiguities for polypeptides.
* Drude patch DISU is not included due to a lack of support in ParmEd and
  OpenMM.
* Thulium(III) ion from `toppar/stream/misc/toppar_ions_won.str` (TM3P) is
  excluded because it collides with the identically-named 4'-methyl,3'-phosphate
  tetrahydrofuran residue.

If you need parameters from these files that were excluded, you may download the
CHARMM TOPPAR files, manually edit the offending entries, and regenerate the
OpenMM XML files at your own risk.

Since the full force field `charmm36_nowaters.xml` is somewhat large, you can
load the `charmm36_protein.xml` force field for protein parameters only, if
applicable.

Due to some other issues, some residues and patches have been split out into
separate files.  Residue and patch collisions may result from using these files,
and you may need to use the `residueTemplates` option of OpenMM's
[ForceField.createSystem()](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField.createSystem)
to successfully assign residues when using these files.  They may be loaded
after the primary force field files:

* For `charmm36_nowaters.xml`: D-amino acids and other modified amino acids are
  in `charmm36_d_modified.xml`, extra carbohydrate patches are in
  `charmm36_carb_imlab.xml`, and CGenFF small molecules are in
  `charmm36_cgenff.xml`.
* For `charmm36_protein.xml`: D-amino acids are in `charmm36_protein_d.xml`.
* For `charmm_polar_2023.xml` (Drude force field): D-amino acids are in
  `charmm_polar_2023_d.xml`.

## Converting

Retrieve and unpack the CHARMM files
```
wget -O toppar.tgz http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz
tar -xzf toppar.tgz
(cd toppar/drude; tar -xzf drude_toppar_2023.tgz)
```
Convert force fields (may take several minutes)
```
./convert_charmm.sh
```
Test force fields (may take several minutes)
```
./test_charmm.sh
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

### Notice

Due to various outstanding issues in ParmEd, including:

* [ParmEd/ParmEd#1383](https://github.com/ParmEd/ParmEd/issues/1383)
* [ParmEd/ParmEd#1384](https://github.com/ParmEd/ParmEd/issues/1384)
* [ParmEd/ParmEd#1395](https://github.com/ParmEd/ParmEd/issues/1395)
* [ParmEd/ParmEd#1396](https://github.com/ParmEd/ParmEd/issues/1396)
* [ParmEd/ParmEd#1397](https://github.com/ParmEd/ParmEd/issues/1397)

it is not currently possible to use a prepackaged version of ParmEd or a version
from its GitHub repository to perform the CHARMM conversion.  Until these issues
are fixed, a custom patched version of ParmEd containing workarounds for these
problems specific to this conversion workflow can be found at
[`patched-for-charmm-conversion` branch of the epretti/ParmEd repository](https://github.com/epretti/ParmEd/tree/patched-for-charmm-conversion).

CHARMM docker image
===================

Prerequisites
-------------
* [Docker](https://www.docker.com/products/docker-desktop)
* [CHARMM (free for nonprofit/academic users)](https://brooks.chem.lsa.umich.edu/register/): downloaded as `charmm.tar.gz`

Building the docker image
-------------------------
After starting the Docker daemon, run, *e.g.*,
```
docker build -t test_charmm .
```
(you can replace `test_charmm` with another name) and a Ubuntu-based Docker
image containing the CHARMM binary will be built.  Note that to keep the image
size small, only the binary and needed Fortran runtime libraries will be present
in the final image, not development tools and other files.  Edit the
`Dockerfile` before building to change this.

Running CHARMM
--------------
The CHARMM executable is `/usr/local/bin/charmm` inside the image.

To manually start the docker image (for testing purposes):
```
docker run -i -t test_charmm
```

To run the tests using a dockerized CHARMM, pass
`--charmm-docker-image test_charmm` to `test_charmm.py` or `test_charmm.sh`.
