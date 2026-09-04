[![CI](https://github.com/openmm/openmmforcefields/actions/workflows/ci.yaml/badge.svg)](https://github.com/openmm/openmmforcefields/actions/workflows/ci.yaml)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/openmmforcefields.svg)](https://anaconda.org/conda-forge/openmmforcefields)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/openmmforcefields.svg)](https://anaconda.org/conda-forge/openmmforcefields)
[![DOI](https://zenodo.org/badge/70107487.svg)](https://zenodo.org/badge/latestdoi/70107487)

# OpenMMForceFields: AMBER, CHARMM, OpenFF, and Espaloma force fields for OpenMM

This repository provides support for AMBER, CHARMM, OpenFF, and Espaloma force fields and small molecule parameterization with GAFF, Espaloma, and Open Force Field Toolkit for OpenMM.

## Supported force fields

**AMBER:** All major AMBER force fields distributed with [AmberTools](https://ambermd.org/AmberTools.php) 24.8 from [conda-forge](https://anaconda.org/conda-forge/ambertools/files), as well as all released [GAFF small molecule force fields](http://ambermd.org/antechamber/gaff.html) through 1.81 (GAFF 1.x) and 2.2.20 (GAFF 2.x).

**CHARMM:** Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [July 2024 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml). *Note that this conversion has not yet been fully validated.*

**Open Force Field Initiative force fields:** All distributed [Open Force Field Initiative](http://openforcefield.org) [force fields](https://openforcefield.org/force-fields/force-fields/), including the [`openff-2.x.y` ("Sage")](https://openforcefield.org/force-fields/force-fields/) and [`smirnoff99Frosst`](https://github.com/openforcefield/smirnoff99Frosst/) series of force fields available through the [`openff-forcefields`](http://github.com/openforcefield/openff-forcefields) package.

**Espaloma:** Currently [`espaloma-0.3.2`](https://github.com/choderalab/espaloma/releases/tag/0.3.2) is supported. See our [first espaloma paper](https://arxiv.org/abs/2010.01196) and our second paper which focuses on [protein-ligand systems and beyond](https://arxiv.org/abs/2307.07085).

## Installation

The OpenMMForceFields package provides additional AMBER and CHARMM biopolymer force fields, small molecule support through GAFF and the [Open Force Field toolkit](http://openforcefield.org), and force field conversion tools.

The easiest way to install this package and its requisite dependencies is via [`conda`](https://conda.io):
```bash
conda install --yes -c conda-forge openmmforcefields
```

## Documentation

Documentation for the AMBER and CHARMM force fields, as well as for the
OpenMMForceFields API for using SMIRNOFF, GAFF, and Espaloma force fields, is
available at <https://openmm.github.io/openmmforcefields>.

### Converting AMBER and CHARMM to OpenMM ffxml files

See the corresponding directories for information on how to use the provided conversion tools:

* `amber/` - AMBER force fields and conversion tools
* `charmm/` - CHARMM force fields and conversion tools

## Frequently Asked Questions (FAQ)

**Q:** What is the minimum version of OpenMM required to use this package?
<br>
**A:** You need at least OpenMM 8.5.1 to use the OpenMMForceFields package.

**Q:** Do you support the new [Amber ff19SB protein force field](https://doi.org/10.1021/acs.jctc.9b00591)?
<br>
**A:** ff19SB and phosaa19SB have now been added to OpenMMForceFields.

**Q:** Do you plan to support other small molecule force fields?
<br>
**A:** If there are other free and open source conda-installable tools for generating parameters for other AMBER- or CHARMM-compatible force fields, we are happy to add support for them!

# [Changelog](https://github.com/openmm/openmmforcefields/releases)

## 0.16.0 Template generator improvements

This release adds support to `SMIRNOFFTemplateGenerator` for virtual sites and constraints in SMIRNOFF force fields, as well as for loading multiple SMIRNOFF force field files into one template generator.  It also adds support for parameterizing molecules spanning more than one residue in an OpenMM `Topology`.  In addition, this release improves the performance of residue template matching and caching.

Note that this release changes the behavior of template generators when no `forcefield` argument is provided.  Previously, the latest supported force field for a given template generator, which could change from release to release, would be selected automatically.  Starting with this release of OpenMMForceFields, specifying the `forcefield` argument explicitly is mandatory for all template generators.

With this release, OpenMMForceFields now requires an OpenMM version no earlier than 8.5.1.

## 0.15.1 More Force Field Updates!

We added Lipid21 in [PR #390](https://github.com/openmm/openmmforcefields/pull/390) and updated the default SMIRNOFF force field to be `openff-2.2.1` in [PR #417](https://github.com/openmm/openmmforcefields/pull/417) 
We removed the ff19ipq force field since it is not supported by Amber, see [PR #419](https://github.com/openmm/openmmforcefields/pull/419) for details.
We now support any barostat in SystemGenerator, see [PR #414](https://github.com/openmm/openmmforcefields/pull/414) for details.

## 0.15.0 Force Field Updates and More!

All template generators now take `template_generator_kwargs` as an optional init parameter.
This removes the `residue_atoms` argument from the `generate_residue_template` method, see [PR #391](https://github.com/openmm/openmmforcefields/pull/391) for details.
This release updates the CHARMM force field to the July 2024 release.
We also updated Amber force fields with latest versions from AmberTools 24.
Compatibility with esploma force fields has also been improved.

For more details, see our [0.15.0 release page](https://github.com/openmm/openmmforcefields/releases/tag/0.15.0) for more details.

## 0.14.1 Bring back GAFFTemplateGenerator for OpenMM >=7.6.0

This release brings back GAFF force feild support for all versions of OpenMM previously supported.
Additionally, we now use the output of parmchk2 for all GAFF parameters.
Previously we used gaff.dat + parmchk2 output to generate forcefield parameters.
Functionally this doesn't change the end user experience but means we do not need to create new forcefield XML files for newer GAFF versions and now support whatever GAFF versions that parmchk2 supports for the installed AmberTools version.
The XML files in `openmmforcefields/ffxml/amber/gaff/ffxml` may be removed in a future release.

## 0.14.0 Bring back GAFFTemplateGenerator

This release effectively reverts the changes in 0.13.0. This release is **only** compatible with OpenMM 8.1.2. No other changes were made.

## 0.13.0 Temporarily remove GAFFTemplateGenerator

This release temporarily removes GAFFTemplateGenerator because of packaging incompatibilities with
AmberTools 23. This functionality is planned to be re-introduced in 0.14.0.

This release is expected to work with Python 3.10-3.12.

Other changes include
* The default force field of `SystemGenerator` was updated from `openff-1.0.0` (code name Parsley) to
  `openff-2.0.0` (code name Sage).

## 0.12.0 Updates for espaloma and support a offxml string in SystemGenerator
See our [0.12.0 release page](https://github.com/openmm/openmmforcefields/releases/tag/0.12.0) for more details.

## 0.11.0 Support for espaloma small molecule parameters
This release adds support for using [espaloma](https://github.com/choderalab/espaloma) to apply small molecule parameters.

* [(PR #182)](https://github.com/openmm/openmmforcefields/pull/179) Add support for espaloma small molecule parameters

## 0.10.0 Updates for OpenMM 7.6 and AMBER GLYCAM addition
This release adds support for the AMBER GLYCAM force field supporting glycans and updates imports for OpenMM 7.6.

* [(PR #156)](https://github.com/openmm/openmmforcefields/pull/156) Add support for GLYCAM
* [(PR #173)](https://github.com/openmm/openmmforcefields/pull/173) Update to OpenMM 7.6 imports

## 0.9.0 Updates for openforcefield 0.9.0 toolkit
This release utilizes the new [openforcefield 0.9.0 toolkit](https://open-forcefield-toolkit.readthedocs.io/en/0.9.0/) now distributed through [conda-forge](https://conda-forge.org/).

This release contains updated CHARMM and AMBER force fields for use with [OpenMM 7.5.0](https://github.com/openmm/openmm/releases/tag/7.5.0) and the new [openforcefield 0.9.0 toolkit](https://open-forcefield-toolkit.readthedocs.io/en/0.9.0/), both now distributed through [conda-forge](https://conda-forge.org/).

* Amber force fields were updated to versions distributed with [AmberTools 20.15](https://anaconda.org/conda-forge/ambertools/files)
* Added AMBER `phosaa14SB` parameters for phosphorylated amino acids
* CHARMM force fields were updated to [July 2020 CHARMM additive force field release](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm)

## 0.8.0 Updates for openforcefield 0.7.1 toolkit
* [(PR #128)](https://github.com/openmm/openmmforcefields/pull/128) Update README for openff-1.2.0 and use openforcefield 0.7.1 toolkit API for identifying installed force fields

## 0.7.5 Bugfix release
* [(PR #127)](https://github.com/openmm/openmmforcefields/pull/127) Fixes a bug where the wrong path was imported for logging; improves docstrings.

## 0.7.4 Bugfix release to ensure compatibility with openforcefield toolkit 0.7.0
* [(PR #121)](https://github.com/openmm/openmmforcefields/pull/121) Add compatibility with [`openforcefield 0.7.0`](https://github.com/openforcefield/openff-toolkit/releases/tag/0.7.0)

## 0.7.3 Bugfix release: Compatibility with openforcefield toolkit 0.7.0 and auto-detection of installed openforcefield force fields
* [(PR #119)](https://github.com/openmm/openmmforcefields/pull/119) Handle `None` partial charges in openforcefield `Molecule` objects (needed in `openforcefield` toolkit 0.7.0)
* [(PR #120)](https://github.com/openmm/openmmforcefields/pull/120) Auto-detect installed SMIRNOFF force fields

## 0.7.2 Bugfix release: More error checking; OpenMM 7.4.2 minimum version requirement

* Raise a `ValueError` if `SystemGenerator` receives a `nonbondedMethod` key in `forcefield_kwargs`; these should go into `periodic_forcefield_kwargs` or `nonperiodic_forcefield_kwargs`.

## 0.7.1 Bugfix release: Fix GAFF AM1-BCC charging bug for some molecules

* When using the OpenEye toolkit, some molecules failed to charge with GAFF. See https://github.com/openforcefield/openff-toolkit/issues/492
* Removed most `perses_jacs_systems`, updating the remaining ones with inputs used in AMBER-TI publication, https://pubs.acs.org/doi/10.1021/acs.jcim.9b00105
* Fix a bug where some molecules with pyramidal atoms would cause exceptions when read from cache

## 0.7.0 User-specified partial charges, SystemGenerator support for periodic and non-periodic topologies, and minor bugfixes

* If `Molecule` objects contain nonzero partial charges, these are used instead of generating new partial charges
* Fix bug in default `periodic_forcefield_kwargs` and `nonperiodic_forcefield_kwargs`

## 0.6.1 Updated README and minor bugfixes

* Fix examples in the `README.md`
* Fix `GAFFTemplateGenerator.gaff_major_version`
* Fix incorrect default SMIRNOFF force field, which is now `openff-1.0.0` (was previously `smirnoff99Frosst-1.1.0`)
* Add `SystemGenerator.SMALL_MOLECULE_FORCEFIELDS` convenience property to list available small molecule force fields
* `SystemGenerator` API changed to support both periodic and non-periodic `Topology` objects for the same generator

## 0.6.0 Updated AMBER force fields (AmberTools 19.9) and small molecule support via GAFF

This release provides updated support for AMBER biopolymer force fields (from AmberTools 19.9) and small molecule support with GAFF 1.x and 2.x, along with experimental support for the new Open Force Field Initiative SMIRNOFF force fields.

**AMBER:** All major AMBER force fields distributed with [AmberTools](https://ambermd.org/AmberTools.php) 19.9 (except ff19SB---see FAQ below), as well as all released [GAFF small molecule force fields](http://ambermd.org/antechamber/gaff.html) through 1.81 (GAFF 1.x) and 2.11 (GAFF 2.x).

**CHARMM:** Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [Aug 2015 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml). *Note that this conversion has not yet been fully validated.*

**Open Force Field Initiative force fields:** All distributed [Open Force Field Initiative](http://openforcefield.org) force fields, including the `smirnoff99Frosst` series and [`openff-1.0.0` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/). *This support is experimental since it requires a development version of OpenMM 7.5.0.*

Residue template generators are provided for both GAFF (`GAFFTemplateGenerator`) and SMIRNOFF (`SMIRNOFFTemplateGenerator`).

This release also contains an experimental new `SystemGenerator` for managing biopolymer and small molecule force field `System` creation via a unified API.

## 0.5.0 Force fields for OpenMM 7.3.1

This release contains updated CHARMM and AMBER force fields for use with OpenMM 7.3.1.

* Amber force fields were updated to versions distributed with [AmberTools 18.0](https://anaconda.org/omnia/ambertools/files)
* Release version metadata in Amber force fields was corrected
* CHARMM force fields were updated to [July 2018 CHARMM additive force field release](http://mackerell.umaryland.edu/charmm_ff.shtml#charmm)
* Experimental Amber GAFF residue template generator released (requires [OpenEye Toolkit](https://docs.eyesopen.com/toolkits/python/index.html))

## 0.4.0 Force fields for OpenMM 7.3.0

This release contains updated CHARMM and AMBER force fields distributed with OpenMM 7.3.0.

Amber force fields were converted from the AmberTools 18 package, while CHARMM force fields were converted from the July 2016 update.

## 0.3.0 Force fields for OpenMM 7.2.0 rc2

This release contains the force fields distributed with OpenMM 7.2.0 rc2

## 0.2.0 Forcefields for OpenMM 7.2.0 rc1

This release contains the force fields distributed with OpenMM 7.2.0 rc1

## 0.1.0 Conda-installable AMBER force fields

This prerelease allows installation of AmberTools 16 via conda.
