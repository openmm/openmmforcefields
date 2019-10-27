[![Build Status](https://travis-ci.org/choderalab/openmm-forcefields.svg?branch=master)](https://travis-ci.org/choderalab/openmm-forcefields?branch=master)
[![DOI](https://zenodo.org/badge/70107487.svg)](https://zenodo.org/badge/latestdoi/70107487)

# CHARMM and AMBER force fields for OpenMM

This repository contains additional CHARMM and AMBER force fields---along with their conversion tools---for OpenMM.

Force field conversion tools make use of [ParmEd](https://github.com/parmed/parmed).

# Using the force fields

## Installation

You can install the additional force fields and conversion tools via [`conda`](https://conda.io):
```bash
conda install --yes -c conda-forge -c omnia openmm-forcefields
```

## Using the CHARMM and AMBER force fields

Once installed, the AMBER force fields will be registered in the `amber/` relative path searched by [`simtk.openmm.app.ForceField`](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).
For example, to specifically access the `ff99SBildn` force field, prepend the `amber` prefix and the `.xml` suffix:
```python
forcefield = ForceField('amber/ff99SBildn.xml')
```
Similarly, the CHARMM force fields will be registered in the `charmm/` relative path.
For example, model system small molecule templates corresponding to amino acids can be accessed with:
```python
forcefield = ForceField('charmm/toppar_all36_prot_model.xml')
```

## Using the AMBER GAFF residue template generator for small molecules

A [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for `ForceField` that automatically generates residue templates for small molecules lacking them using [GAFF or GAFF2](http://ambermd.org/antechamber/gaff.html) is also provided.
This version supports the [OpenEye Toolkit](https://docs.eyesopen.com/toolkits/python/index.html) as a means to specify the small molecules to be matched to `Residue`s in the OpenMM `Topology` object.

For example, to use GAFF, first create a `ForceField` object using `gaff.xml`, and then add a residue template generator created with [`OEMol`](https://docs.eyesopen.com/toolkits/python/oechemtk/OEChemClasses/OEMol.html) versions of your molecules of interest:
```python
# Create a ForceField object with the GAFF parameters
forcefield = ForceField('amber/gaff.xml')
# Create a residue template generator factory, and use it to create a template generator that knows about the specified OEMol molecules (oemols)
# We ask for lazy parameterization (parameters are generated only when needed) and the latest GAFF version 1
from openmmforcefields.gaff import TemplateGeneratorFactory
template_generator = TemplateGeneratorFactory.create_generator(oemols, lazy=True, gaff_version='1')
# Register the new template generator
forcefield.registerTemplateGenerator(template_generator)
```
By default, parameters are generated lazily.
If we want all molecules to be parameterized and registered immediately (which may be slow), we can set `lazy=False`:
```python
template_generator = TemplateGeneratorFactory.create_generator(oemols, lazy=False, gaff_version='1')
```
Note that only small molecules that correspond to single residues are supported at this time---this does not yet work for biopolymer residues.

We can also select GAFF 2 if we desire:
```python
template_generator = TemplateGeneratorFactory.create_generator(oemols, gaff_version='2')
```

The [canonical AM1-BCC charging method](https://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html) is used to assign ELF10 charges.

# Converting CHARMM and AMBER force fields

See the corresponding directories for information on how to use the provided conversion tools:

* `charmm/` - CHARMM force fields and conversion tools
* `amber/` - AMBER force fields and conversion tools

# Changelog

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
