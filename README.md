[![Build Status](https://travis-ci.org/choderalab/openmm-forcefields.svg?branch=master)](https://travis-ci.org/choderalab/openmm-forcefields?branch=master)
[![DOI](https://zenodo.org/badge/70107487.svg)](https://zenodo.org/badge/latestdoi/70107487)

# AMBER and CHARMM force fields for OpenMM

This repository adds support for small molecule parameterization and additional AMBER and CHARMM force fields for OpenMM.

## Supported force fields

*AMBER:* All major AMBER force fields distributed with AmberTools 18.0, as well as the GAFF small molecule force fields through 1.81 and 2.11
*CHARMM:* Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [Aug 2015 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml)
*Open Force Field:* All distributed [Open Force Field Initiative](http://openforcefield.org) force fields through `openforcefield-1.0.0` ("Parsley"). Newer parameter sets can be obtained by updating the [`openforcefields` package](https://github.com/openforcefield/openforcefields) with `conda`.

# Using the force fields

## Installation

You can install the additional AMBER and CHARMM force fields and conversion tools provided by `openmm-forcefields` via [`conda`](https://conda.io):
```bash
conda install --yes -c conda-forge -c omnia openmm-forcefields
```
This will also install all the dependencies needed to generate small molecule parameters using the [`openforcefield` toolkit](http://openforcefield.org).

## Using the AMBER and CHARMM force fields

Once installed, the AMBER force fields will be registered in the `amber/` relative path searched by [`simtk.openmm.app.ForceField`](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).
For example, to specifically access the new `ff19SB` force field, prepend the `amber` prefix and the `.xml` suffix:
```python
forcefield = ForceField('amber/ff19SB.xml')
```
or to specifically access the `ff99SBildn` force field:
```python
forcefield = ForceField('amber/ff99SBildn.xml')
```

Similarly, the CHARMM force fields will be registered in the `charmm/` relative path.
For example, model system small molecule templates corresponding to amino acids can be accessed with:
```python
forcefield = ForceField('charmm/toppar_all36_prot_model.xml')
```

## Using AMBER GAFF 1/2 for small molecules

The `openmm-forcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using [GAFF](http://ambermd.org/antechamber/gaff.html) versions 1 or 2.

The [`openforcefield` toolkit](http://openforcefield.org) is used to provide an interface with cheminformatics toolkits to interact with [`antechamber`](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) package to generate parameters for small molecules.
By default, the [`openforcefield` toolkit](http://github.com/openforcefield/openforcefield) will make use of the free and open source [RDKit cheminformatics toolkit](https://www.rdkit.org/) that is installed automatically, but will optionally use the [OpenEye toolkit](https://docs.eyesopen.com/toolkits/python/index.html) if it is installed and licensed.
The OpenEye toolkit is available [for free for academics for non-IP generating academic research](https://www.eyesopen.com/academic-licensing).

Generation of OpenMM-compatible parameters is handled through `openmmforcefields.generators.GAFFTemplateGenerator`.
Because the [OpenMM `Topology` object](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.topology.Topology.html#simtk.openmm.app.topology.Topology) used by [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) does not know the precise chemical identity of molecules represented in the topology---which contain only elements and bonds between them, without stereochemical or bond order information---it is necessary to instruct `GAFFTemplateGenerator` which small molecules it will encounter in processing the `Topology` object ahead of time; it then matches these by element and bond pattern.

To do this, it is necessary to specify one or more [`openforcefield.topology.Molecule`](https://open-forcefield-toolkit.readthedocs.io/en/0.5.1/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) objects which can easily be created from many different representations of small molecules, including SMILES strings and common molecule storage formats.
There are many ways to create an [openforcefield `Molecule` object](https://open-forcefield-toolkit.readthedocs.io/en/0.5.1/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) from various file formats as well---see the [API docs](https://open-forcefield-toolkit.readthedocs.io/en/0.5.1/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) for more details.

The `openforcefield` toolkit charging method [`Molecule.compute_partial_charges_am1bcc`](https://open-forcefield-toolkit.readthedocs.io/en/0.5.1/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule.compute_partial_charges_am1bcc) is used to assign partial charges.
The [canonical AM1-BCC charging method](https://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html) is used to assign ELF10 charges if the OpenEye Toolkit is available.
If not, [`antechamber -c bcc`](http://ambermd.org/antechamber/) from the [AmberTools]](http://ambermd.org/AmberTools.php) distribution (which uses the `sqm` semiempirical quantum chemical package) is used to assign AM1-BCC charges.

*Note:* The molecule you specify must have the all protons and stereochemistry explicitly specified, and must match the exact protonation and tautomeric state of the molecule that will be found in your OpenMM `Topology` object.

`GAFFTemplateGenerator` also supports the ability to specify a JSON cache file, allowing parameters for small molecules to be computed only once and then cached thereafter.

### Examples using `GAFFTemplateGenerator` to generate small molecule GAFF parameters

Create a template generator for GAFF for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
from openforcefield.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
from openmoltools.forcefield_generators import GAFFTemplateGenerator
gaff = GAFFTemplateGenerator(molecules=molecule, gaff_version='1.81')
from simtk.openmm.app import ForceField
forcefield = ForceField(gaff.gaff_xml_filename, 'amber14-all.xml', 'tip3p.xml')
forcefield.registerTemplateGenerator(gaff)
```

Create a template generator for GAFF2 for multiple molecules read from a file:
```python
molecules = Molecule.from_file('molecules.sdf')
gaff = OEGAFFTemplateGenerator(molecules=molecules, gaff_version='2.11')
```
You can also add molecules later on after the generator has been registered:
```python
gaff.add_molecules(molecule)
gaff.add_molecules([molecule1, molecule2])
```
To check which GAFF versions are supported, check the `SUPPORTED_GAFF_VERSIONS` attribute:
```python
>>> print(GAFFTemplateGenerator.SUPPORTED_GAFF_VERSIONS)
['1.4', '1.8', '1.81', '2.1', '2.11']
```
You can optionally create or use a tiny database cache of pre-parameterized molecules:
```python
gaff = GAFFTemplateGenerator(cache='gaff-molecules.json', gaff_version='1.80')
```
Newly parameterized molecules will be written to the cache, saving time next time!

## Automating force field management with `SystemGenerator`s

The `openmm-forcefields` package provides several `openmmforcefields.generators.SystemGenerator` subclasses that handle management of common force fields transparently for you.

### Examples using `GAFFSystemGenerator` to automate the use of AMBER force fields with GAFF for small molecule parameterization

Here's an example that uses GAFF 2.11 along with the new `ff19SB` generation of AMBER force fields (and compatible solvent models):
```python
from openmmforcefields.generators import GAFFSystemGenerator
system_generator = GAFFSystemGenerator(amber_forcefield='ff19SB', solvent_model='tip3p', gaff_version='2.11', molecules=molecules, cache='gaf-2.11.json')
```

# Frequently Asked Questions (FAQ)

**Q:** What is the minimum version of OpenMM required?

**A:** OpenMM 7.4.1 is the minimum required version.

# Converting AMBER and CHARMM force fields to OpenMM ffxml files

See the corresponding directories for information on how to use the provided conversion tools:

* `amber/` - AMBER force fields and conversion tools
* `charmm/` - CHARMM force fields and conversion tools

# Changelog

## 0.6.0 Force fields for OpenMM 7.4.1 and GAFF support

This release contains updated CHARMM and AMBER force fields for use with OpenMM 7.4.1.

* Amber force fields were updated to the versions distributed with [AmberTools 19.9](http://ambermd.org/AmberTools.php), which includes the [Amber ff19SB protein force field](https://chemrxiv.org/articles/ff19SB_Amino-Acid_Specific_Protein_Backbone_Parameters_Trained_Against_Quantum_Mechanics_Energy_Surfaces_in_Solution/8279681/1).
* Support for GAFF via [`antechamber`](http://ambermd.org/antechamber/) and the [`openforcefield` toolkit](http://openforcefield.org).

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
