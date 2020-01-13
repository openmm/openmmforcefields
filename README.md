[![Build Status](https://travis-ci.org/openmm/openmm-forcefields.svg?branch=master)](https://travis-ci.org/openmm/openmm-forcefields?branch=master)
[![DOI](https://zenodo.org/badge/70107487.svg)](https://zenodo.org/badge/latestdoi/70107487)

# AMBER and CHARMM force fields for OpenMM

This repository provides support for AMBER and CHARMM force fields and small molecule parameterization with GAFF and the Open Force Field Toolkit for OpenMM.

## Supported force fields

**AMBER:** All major AMBER force fields distributed with [AmberTools](https://ambermd.org/AmberTools.php) 19.9 (except ff19SB---see FAQ below), as well as all released [GAFF small molecule force fields](http://ambermd.org/antechamber/gaff.html) through 1.81 (GAFF 1.x) and 2.11 (GAFF 2.x).

**CHARMM:** Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [Aug 2015 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml). *Note that this conversion has not yet been fully validated.*

**Open Force Field Initiative force fields:** All distributed [Open Force Field Initiative](http://openforcefield.org) force fields, including the `smirnoff99Frosst` series and [`openff-1.0.0` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/). *This support is experimental since it requires a development version of OpenMM 7.5.0.*

# Using the force fields

## Installation

The `openmm-forcefields` package provides additional AMBER and CHARMM biopolymer force fields, small molecule support through GAFF and the [Open Force Field toolkit](http://openforcefield.org), and force field conversion tools.

### Hacky temporary installation workaround

```bash
# Install dependencies (assuming you have conda-forge, omnia, and optionally openeye channels already added)
conda install --yes pyyaml "ambertools>=18.0" lxml networkx tinydb "openforcefield>=0.6.0" "openforcefields>=1.0.0"
# Install openmm-forcefields
pip install git+https://github.com/openmm/openmm-forcefields.git
# Patch your OpenMM installation to be able to use SMIRNOFF (assumes OpenMM 7.4.1)
wget -q https://raw.githubusercontent.com/openmm/openmm/4f48402f1be3e0f049ae0e8595db638d297b0d75/wrappers/python/simtk/openmm/app/forcefield.py -O `python -c "from simtk.openmm.app import forcefield; print(forcefield.__file__)"`
```

### NOTE: The method below does not work yet---we are working on building this conda package

To install this package, you will need to

The easiest way to install this package and its requisite dependencies is via [`conda`](https://conda.io):
```bash
conda install --yes -c conda-forge -c omnia openmm-forcefields
```

## Using the AMBER and CHARMM biopolymer force fields

This repository contains force fields for use with the [OpenMM `ForceField` class](http://docs.openmm.org/latest/userguide/application.html#force-fields) for parameterizing biomolecular systems.
If you're not familiar with this approach to applying parameters to biomolecular systems, please see the [OpenMM User Guide](http://docs.openmm.org/latest/userguide/application.html#force-fields).

### Using the AMBER force fields

Once installed, the AMBER force fields will be registered in the `amber/` relative path searched by [`simtk.openmm.app.ForceField`](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField).

For example, to specify the newer recommended [`ff14SB`](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255) force field and accompanying recommended ions and solvent models (corresponding to force fields loaded in LEaP with `leaprc.protein.ff14SB`), prepend the `amber` prefix and the `.xml` suffix:
```python
forcefield = ForceField('amber/protein.ff14SB.xml')
```
To access just the [`ff14SB`](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255) force field converted from `oldff/leaprc.ff14SB`:
```python
forcefield = ForceField('amber/ff14SB.xml')
```
or to specifically access the older (now outdated) `ff99SBildn` force field converted from `oldff/leaprc.ff14SB`:
```python
forcefield = ForceField('amber/ff99SBildn.xml')
```
The TIP3P conversion also includes the [Joung and Cheatham recommended salt models](https://doi.org/10.1021/jp8001614) (`parm/frcmod.ionsjc_tip3p`) and recommended divalent counterion parameters (`parm/frcmod.ions234lm_126_tip3p`):
```python
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
```

### Using the CHARMM force fields

Similarly, the CHARMM force fields will be registered in the `charmm/` relative path.
For example, model system small molecule templates corresponding to amino acids can be accessed with:
```python
forcefield = ForceField('charmm/toppar_all36_prot_model.xml')
```

## Using AMBER GAFF 1.x and 2.x for small molecules

The `openmm-forcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using [GAFF](http://ambermd.org/antechamber/gaff.html) versions 1 or 2.

### Cheminformatics toolkits

The [`openforcefield` toolkit](http://openforcefield.org) is used to provide an interface with cheminformatics toolkits to interact with [`antechamber`](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) package to generate parameters for small molecules.
By default, the [`openforcefield` toolkit](http://github.com/openforcefield/openforcefield) will make use of the free and open source [RDKit cheminformatics toolkit](https://www.rdkit.org/) that is installed automatically, but will optionally use the [OpenEye toolkit](https://docs.eyesopen.com/toolkits/python/index.html) if it is installed and licensed.
The OpenEye toolkit is available [for free for academics for non-IP-generating academic research](https://www.eyesopen.com/academic-licensing).

### On-the-fly template generation for small molecules

Generation of OpenMM-compatible parameters for small molecules encountered in an OpenMM `Topology` is handled through `openmmforcefields.generators.GAFFTemplateGenerator`.
Because the [OpenMM `Topology` object](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.topology.Topology.html#simtk.openmm.app.topology.Topology) used by [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) does not know the precise chemical identity of molecules represented in the topology---which contain only elements and bonds between them, without stereochemical or bond order information---it is necessary to instruct `GAFFTemplateGenerator` which small molecules it will encounter in processing the `Topology` object ahead of time; it then matches these by element and bond pattern.

### Specifying molecules

To do this, it is necessary to specify one or more [`openforcefield.topology.Molecule`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) objects which can easily be created from many different representations of small molecules, including SMILES strings and common molecule storage formats.
There are many ways to create an [openforcefield `Molecule` object](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) from various file formats as well---see the [API docs](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule) for more details.

The `openforcefield` toolkit charging method [`Molecule.compute_partial_charges_am1bcc`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openforcefield.topology.Molecule.html#openforcefield.topology.Molecule.compute_partial_charges_am1bcc) is used to assign partial charges.
The [canonical AM1-BCC charging method](https://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html) is used to assign ELF10 charges if the OpenEye Toolkit is available.
If not, [Antechamber](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) distribution (which uses the `sqm` semiempirical quantum chemical package) is used to assign AM1-BCC charges (`antechamber -c bcc`).

*Note:* The `Molecule` object must have the all protons and stereochemistry explicitly specified, and must match the exact protonation and tautomeric state of the molecule that will be found in your OpenMM `Topology` object.
The atom ordering need not be the same.

### Caching

`GAFFTemplateGenerator` also supports the ability to specify a cache filename, allowing parameters for small molecules to be computed only once and then cached in the specified cache file thereafter.

### Examples using `GAFFTemplateGenerator` to generate small molecule GAFF parameters

Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an openforcefield Molecule object for benzene from SMILES
from openforcefield.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
# Create the GAFF template generator
from openmmforcefields.generators import GAFFTemplateGenerator
gaff = GAFFTemplateGenerator(molecules=molecule)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from simtk.openmm.app import ForceField
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
# Register the GAFF template generator
forcefield.registerTemplateGenerator(gaff.generator)
# You can now parameterize an OpenMM Topology object that contains the specified molecule.
# forcefield will load the appropriate GAFF parameters when needed, and antechamber
# will be used to generate small molecule parameters on the fly.
from simtk.openmm.app import PDBFile
pdbfile = PDBFile('t4-lysozyme-L99A-with-benzene.pdb')
system = forcefield.createSystem(pdbfile.topology)
```
The latest available GAFF version is used if none is specified.
You can check which GAFF version is in use with
```python
>>> generator.gaff_version
'2.11'
````
Create a template generator for a specific GAFF version for multiple molecules read from an SDF file:
```python
molecules = Molecule.from_file('molecules.sdf')
gaff = GAFFTemplateGenerator(molecules=molecules, forcefield='gaff-2.11')
```
You can also add molecules to the generator later, even after the generator has been registered:
```python
gaff.add_molecules(molecule)
gaff.add_molecules([molecule1, molecule2])
```
To check which GAFF versions are supported, examine the `INSTALLED_FORCEFIELDS` attribute:
```python
>>> print(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)
['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11']
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
generator = GAFFTemplateGenerator(cache='gaff-molecules.json', forcefield='gaff-1.80')
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Using the Open Force Field Initiative SMIRNOFF small molecule force fields

The `openmm-forcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.forcefield.ForceField.html#simtk.openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using the [Open Force Field Initiative](http://openforcefield.org) [SMIRNOFF](https://open-forcefield-toolkit.readthedocs.io/en/0.6.0/smirnoff.html) small molecule force fields.
This includes the [`openff-1.0.0` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/) small molecule force field.

The `SMIRNOFFTemplateGenerator` residue template generator operates in a manner very similar to `GAFFTemplateGenerator`, so we only highlight its differences here.

### Examples using `SMIRNOFFTemplateGenerator` to generate small molecule SMIRNOFF parameters

Create a SMIRNOFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an openforcefield Molecule object for benzene from SMILES
from openforcefield.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
# Create the SMIRNOFF template generator with the default installed force field (openff-1.0.0)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from simtk.openmm.app import ForceField
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)
```
The latest official Open Force Field Initiative release ([`openff-1.0.0` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/)) is used if none is specified.
You can check which SMIRNOFF force field is in use with
```python
>>> smirnoff.smirnoff_filename
'/Users/choderaj/Library/Caches/Python-Eggs/openforcefields-1.0.0-py3.7.egg-tmp/openforcefields/offxml/openff-1.0.0.offxml'
```
Create a template generator for a specific SMIRNOFF force field for multiple molecules read from an SDF file:
```python
molecules = Molecule.from_file('molecules.sdf')
# Create a SMIRNOFF residue template generator from the official openff-1.0.0 release,
# which is installed automatically
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='openff-1.0.0')
# Create a SMIRNOFF residue template generator from the official smirnoff99Frosst-1.1.0 release,
# which is installed automatically
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='smirnoff99Frosst-1.1.0')
# Use a local .offxml file instead
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='local-file.offxml')
```
You can also add molecules to the generator later, even after the generator has been registered:
```python
smirnoff.add_molecules(molecule)
smirnoff.add_molecules([molecule1, molecule2])
```
To check which SMIRNOFF force fields are automatically installed, examine the `INSTALLED_FORCEFIELDS` attribute:
```python
>>> print(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)
['smirnoff99Frosst-1.1.0', 'openff-1.0.0']
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
smirnoff = SMIRNOFFTemplateGenerator(cache='smirnoff-molecules.json', forcefield='openff-1.0.0')
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Automating force field management with `SystemGenerator`

The `openmm-forcefields` package provides the `openmmforcefields.generators.SystemGenerator` class that handles management of common force fields transparently for you.

### Using `SystemGenerator` to automate the use of AMBER force fields with GAFF for small molecule parameterization

Here's an example that uses GAFF 2.11 along with the new `ff14SB` generation of AMBER force fields (and compatible solvent models) to generate an OpenMM `System` object from an [Open Force Field `Topology`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openforcefield.topology.Topology.html#openforcefield.topology.Topology) object:
```python
# Define the keyword arguments to feed to ForceField
from simtk import unit
from simtk.openmm import app
forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
# Initialize a SystemGenerator using GAFF
from openmmforcefields.generators import SystemGenerator
system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'], small_molecule_forcefield='gaff-2.11', forcefield_kwargs=forcefield_kwargs, cache='db.json')
# Create an OpenMM System from an Open Force Field toolkit Topology object
system = system_generator.create_system(openforcefield_topology)
# Alternatively, create an OpenMM System from an OpenMM Topology object and a list of openforcefield Molecule objects
system = system_generator.create_system(openmm_topology, molecules=molecules)
```
Parameterized molecules are cached in `db.json`.
Parameters for multiple force fields can be held in the same cache file.

By default, `SystemGenerator` will use `PME` for periodic systems and `NoCutoff` for non-periodic systems.
You can modify this behavior with the optional `periodic_nonbonded_method` and `nonperiodic_nonbonded_method` arguments:
```python
from simtk.openmm import app
system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'], periodic_nonbonded_method=app.LJPME, nonperiodic_nonbonded_method=app.CutoffNonPeriodic)
```

To use the [Open Force Field `openff-1.0.0` ("Parsley") force field](https://openforcefield.org/news/introducing-openforcefield-1.0/) instead of GAFF 2.11, we would have instead specified `small_molecule_forcefield='openff-1.0.0'`.

# Frequently Asked Questions (FAQ)

**Q:** What is the minimum version of OpenMM required to use this package?
<br>
**A:** You need at least OpenMM 7.4.1 to use the `openmm-forcefields` package.

**Q:** Do you support the new [Amber ff19SB protein force field](https://chemrxiv.org/articles/ff19SB_Amino-Acid_Specific_Protein_Backbone_Parameters_Trained_Against_Quantum_Mechanics_Energy_Surfaces_in_Solution/8279681/1)?
<br>
**A:** [ParmEd](http://github.com/parmed/parmed), which is used to convert these force fields to OpenMM format, [does not currently support the conversion of AMBER CMAP forces to OpenMM](https://github.com/ParmEd/ParmEd/issues/1066), so we do not yet support this force field, but hope to add support soon.

**Q:** Do you plan to support other small molecule force fields?
<br>
**A:** If there are other free and open source conda-installable tools for generating parameters for other AMBER- or CHARMM-compatible force fields, we are happy to add support for them!

# Converting AMBER and CHARMM to OpenMM ffxml files

See the corresponding directories for information on how to use the provided conversion tools:

* `amber/` - AMBER force fields and conversion tools
* `charmm/` - CHARMM force fields and conversion tools

# Changelog

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
