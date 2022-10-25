[![CI](https://github.com/openmm/openmmforcefields/workflows/CI/badge.svg?branch=master)](https://github.com/openmm/openmmforcefields/actions?query=workflow%3ACI)
[![conda-forge badge](https://anaconda.org/conda-forge/openmmforcefields/badges/installer/conda.svg)](https://anaconda.org/conda-forge/openmmforcefields)
[![conda-forge downloads](https://anaconda.org/conda-forge/openmmforcefields/badges/downloads.svg)](https://anaconda.org/conda-forge/openmmforcefields)
[![DOI](https://zenodo.org/badge/70107487.svg)](https://zenodo.org/badge/latestdoi/70107487)

# AMBER and CHARMM force fields for OpenMM

This repository provides support for AMBER and CHARMM force fields and small molecule parameterization with GAFF and the Open Force Field Toolkit for OpenMM.

## Supported force fields

**AMBER:** All major AMBER force fields distributed with [AmberTools](https://ambermd.org/AmberTools.php) 20.15 from [conda-forge](https://anaconda.org/conda-forge/ambertools/files) (except ff19SB---see FAQ below), as well as all released [GAFF small molecule force fields](http://ambermd.org/antechamber/gaff.html) through 1.81 (GAFF 1.x) and 2.11 (GAFF 2.x).

**CHARMM:** Non-polarizable protein, nucleic acid, and pre-parameterized small molecule force fields available in in the [July 2020 CHARMM36 force field release from the Mackerell website](http://mackerell.umaryland.edu/charmm_ff.shtml). *Note that this conversion has not yet been fully validated.*

**Open Force Field Initiative force fields:** All distributed [Open Force Field Initiative](http://openforcefield.org) [force fields](https://openforcefield.org/force-fields/force-fields/), including the [`openff-1.x.y` ("Parsley")](https://openforcefield.org/force-fields/force-fields/) and [`smirnoff99Frosst`](https://github.com/openforcefield/smirnoff99Frosst/) series of force fields available through the [`openff-forcefields`](http://github.com/openforcefield/openff-forcefields) package. This is now supported in OpenMM 7.5.0 and later.

# Using the force fields

## Installation

The `openmmforcefields` package provides additional AMBER and CHARMM biopolymer force fields, small molecule support through GAFF and the [Open Force Field toolkit](http://openforcefield.org), and force field conversion tools.

The easiest way to install this package and its requisite dependencies is via [`conda`](https://conda.io):
```bash
conda install --yes -c conda-forge openmmforcefields
```

If you optionally have the [OpenEye Toolkits](https://www.eyesopen.com/toolkit-development) installed, `openmmforcefields` will use these to accelerate small molecule parameterization.
Free academic licenses are available for [bona fide academic research](https://www.eyesopen.com/academic-licensing), while licenses for IP generation are available [for a fee](https://www.eyesopen.com/pricing).

Support for the SMIRNOFF residue template or system generators requires OpenMM 7.4.2 or later.

## Using the AMBER and CHARMM biopolymer force fields

This repository contains force fields for use with the [OpenMM `ForceField` class](http://docs.openmm.org/latest/userguide/application.html#force-fields) for parameterizing biomolecular systems.
If you're not familiar with this approach to applying parameters to biomolecular systems, please see the [OpenMM User Guide](http://docs.openmm.org/latest/userguide/application.html#force-fields).

### Using the AMBER force fields

Once installed, the AMBER force fields will be registered in the `amber/` relative path searched by [`openmm.app.ForceField`](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField).

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

The `openmmforcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using [GAFF](http://ambermd.org/antechamber/gaff.html) versions 1 or 2.

### Cheminformatics toolkits

The [`openff-toolkit`](http://openforcefield.org) is used to provide an interface with cheminformatics toolkits to interact with [`antechamber`](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) package to generate parameters for small molecules.
By default, the [`openff-toolkit`](http://github.com/openforcefield/openff-toolkit) will make use of the free and open source [RDKit cheminformatics toolkit](https://www.rdkit.org/) that is installed automatically, but will optionally use the [OpenEye toolkit](https://docs.eyesopen.com/toolkits/python/index.html) if it is installed and licensed.
The OpenEye toolkit is available [for free for academics for non-IP-generating academic research](https://www.eyesopen.com/academic-licensing).

### On-the-fly template generation for small molecules

Generation of OpenMM-compatible parameters for small molecules encountered in an OpenMM `Topology` is handled through `openmmforcefields.generators.GAFFTemplateGenerator`.
Because the [OpenMM `Topology` object](http://docs.openmm.org/latest/api-python/generated/openmm.app.topology.Topology.html#openmm.app.topology.Topology) used by [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) does not know the precise chemical identity of molecules represented in the topology---which contain only elements and bonds between them, without stereochemical or bond order information---it is necessary to instruct `GAFFTemplateGenerator` which small molecules it will encounter in processing the `Topology` object ahead of time; it then matches these by element and bond pattern.

### Specifying molecules

To do this, it is necessary to specify one or more [`openff.toolkit.topology.Molecule`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Molecule.html#openff.toolkit.topology.Molecule) objects which can easily be created from many different representations of small molecules, including SMILES strings and common molecule storage formats.
There are many ways to create an [OpenFF `Molecule` object](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Molecule.html#openff.toolkit.topology.Molecule) from various file formats as well---see the [API docs](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Molecule.html#openff.toolkit.topology.Molecule) for more details.

### Partial charges for small molecules

If the provided molecule(s) contain nonzero (user-specified) partial charges, they will be used.
If they do _not_ contain partial charges, the `openff-toolkit` charging method [`Molecule.compute_partial_charges_am1bcc`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Molecule.html#openff.toolkit.topology.Molecule.compute_partial_charges_am1bcc) is used to assign partial charges.
The [canonical AM1-BCC charging method](https://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html) is used to assign ELF10 charges if the OpenEye Toolkit is available.
If not, [Antechamber](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) distribution (which uses the `sqm` semiempirical quantum chemical package) is used to assign AM1-BCC charges (`antechamber -c bcc`).

*Note:* The `Molecule` object must have the all protons and stereochemistry explicitly specified, and must match the exact protonation and tautomeric state of the molecule that will be found in your OpenMM `Topology` object.
The atom ordering need not be the same.

*Note:* The first time a `Molecule` is specified, added, or cached, if it lacks partial charges, the automatically generated charges will be cached and reused; if it contains user-specified partial charges, those charges will be used and cached. Adding the molecule again with a different set of charges will have no effect on changing which charges are assigned.

### Caching

`GAFFTemplateGenerator` also supports the ability to specify a cache filename, allowing parameters for small molecules to be computed only once and then cached in the specified cache file thereafter.

### Examples using `GAFFTemplateGenerator` to generate small molecule GAFF parameters

Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
# Create the GAFF template generator
from openmmforcefields.generators import GAFFTemplateGenerator
gaff = GAFFTemplateGenerator(molecules=molecule)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
# Register the GAFF template generator
forcefield.registerTemplateGenerator(gaff.generator)
# You can now parameterize an OpenMM Topology object that contains the specified molecule.
# forcefield will load the appropriate GAFF parameters when needed, and antechamber
# will be used to generate small molecule parameters on the fly.
from openmm.app import PDBFile
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

The `openmmforcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using the [Open Force Field Initiative](http://openforcefield.org) [SMIRNOFF](https://open-forcefield-toolkit.readthedocs.io/en/0.6.0/smirnoff.html) small molecule force fields.
This includes the [`openff-1.0.0` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/) small molecule force field, as well as [newer versions of this force field](https://github.com/openforcefield/openff-forcefields).

The `SMIRNOFFTemplateGenerator` residue template generator operates in a manner very similar to `GAFFTemplateGenerator`, so we only highlight its differences here.

### Examples using `SMIRNOFFTemplateGenerator` to generate small molecule SMIRNOFF parameters

Create a SMIRNOFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
# Create the SMIRNOFF template generator with the default installed force field (openff-1.0.0)
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)
```
The latest official Open Force Field Initiative release ([`openff-1.2.0`](https://github.com/openforcefield/openff-forcefields) of the ["Parsley" small molecule force field](https://openforcefield.org/news/introducing-openforcefield-1.0/)) is used if none is specified.
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
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='openff-1.2.0')
# Create a SMIRNOFF residue template generator from the official smirnoff99Frosst-1.1.0 release,
# which is installed automatically
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield='smirnoff99Frosst-1.2.0')
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
['openff-1.0.1', 'openff-1.1.1', 'openff-1.0.0-RC1', 'openff-1.2.0', 'openff-1.1.0', 'openff-1.0.0', 'openff-1.0.0-RC2', 'smirnoff99Frosst-1.0.2', 'smirnoff99Frosst-1.0.0', 'smirnoff99Frosst-1.1.0', 'smirnoff99Frosst-1.0.4', 'smirnoff99Frosst-1.0.8', 'smirnoff99Frosst-1.0.6', 'smirnoff99Frosst-1.0.3', 'smirnoff99Frosst-1.0.1', 'smirnoff99Frosst-1.0.5', 'smirnoff99Frosst-1.0.9', 'smirnoff99Frosst-1.0.7']
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
smirnoff = SMIRNOFFTemplateGenerator(cache='smirnoff-molecules.json', forcefield='openff-1.2.0')
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.


## Using espaloma to generate small molecule force fields

The `openmmforcefields` package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that can automatically generate OpenMM residue templates for small molecules lacking parameters using [espaloma](https://github.com/choderalab/espaloma) via one of its released force fields, provided `espaloma` and its dependencies are installed.
`espaloma` uses a [graph convolutional model](https://arxiv.org/abs/2010.01196) to generate both valence parameters and fast partial charges.

The `EspalomaTemplateGenerator` residue template generator operates in a manner very similar to `GAFFTemplateGenerator`, so we only highlight its differences here.

**This feature is currently experimental and its API is subject to change.**

### Examples using `EspalomaTemplateGenerator` to generate small molecule parameters with espaloma

Create an espaloma template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit.topology import Molecule
molecule = Molecule.from_smiles('c1ccccc1')
# Create the SMIRNOFF template generator with the released espaloma-0.2.0 force field
from openmmforcefields.generators import EspalomaTemplateGenerator
espaloma = EspalomaTemplateGenerator(molecules=molecule, forcefield='espaloma-0.2.2')
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField
forcefield = ForceField('amber/protein.ff14SB.xml', 'amber/tip3p_standard.xml', 'amber/tip3p_HFE_multivalent.xml')
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(espaloma.generator)
```
Create a template generator for a specific espaloma force field for multiple molecules read from an SDF file:
```python
molecules = Molecule.from_file('molecules.sdf')
# Create an espaloma residue template generator from the espaloma-0.2.0 release,
# retrieving it automatically from GitHub release artifacts
espaloma = EspalomaTemplateGenerator(molecules=molecules, forcefield='espaloma-0.2.2')
# Create an espaloma residue template generator from an espaloma model retrieved from a URL
espaloma = EspalomaTemplateGenerator(molecules=molecules, forcefield='https://github.com/choderalab/espaloma/releases/download/0.2.2/espaloma_0.2.2.pt')
# Create an espaloma residue template generator from an espaloma model stored in a local file
espaloma = EspalomaTemplateGenerator(molecules=molecules, forcefield='/path/to/espaloma_0.2.2.pt')
```
You can also add molecules to the generator later, even after the generator has been registered:
```python
smirnoff.add_molecules(molecule)
smirnoff.add_molecules([molecule1, molecule2])
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
espaloma = EspalomaTemplateGenerator(cache='espaloma-molecules.json', forcefield='espaloma-0.2.2')
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Automating force field management with `SystemGenerator`

The `openmmforcefields` package provides the `openmmforcefields.generators.SystemGenerator` class that handles management of common force fields transparently for you.

### Using `SystemGenerator` to automate the use of AMBER force fields with GAFF, OpenFF, or espaloma for small molecule parameterization

Here's an example that uses GAFF 2.11 along with the new `ff14SB` generation of AMBER force fields (and compatible solvent models) to generate an OpenMM `System` object from an [Open Force Field `Topology`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Topology.html#openff.toolkit.topology.Topology) object:
```python
# Define the keyword arguments to feed to ForceField
from openmm import unit
from openmm import app
forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
# Initialize a SystemGenerator using GAFF
from openmmforcefields.generators import SystemGenerator
system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'], small_molecule_forcefield='gaff-2.11', forcefield_kwargs=forcefield_kwargs, cache='db.json')
# Create an OpenMM System from an OpenMM Topology object
system = system_generator.create_system(openmm_topology)
# Alternatively, create an OpenMM System from an OpenMM Topology object and a list of OpenFF Molecule objects
molecules = Molecule.from_file('molecules.sdf', file_format='sdf')
system = system_generator.create_system(openmm_topology, molecules=molecules)
```
Parameterized molecules are cached in `db.json`.
Parameters for multiple force fields can be held in the same cache file.

By default, `SystemGenerator` will use `PME` for periodic systems and `NoCutoff` for non-periodic systems.
You can modify this behavior with the optional `periodic_forcefield_kwargs` and `nonperiodic_forcefield_kwargs` arguments, which are used to update `forcefield_kwargs` depending on whether the system is periodic or non-periodic:
```python
from openmm import app
system_generator = SystemGenerator(forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    periodic_forcefield_kwargs={'nonbondedMethod' : app.LJPME},
    nonperiodic_forcefield_kwargs={'nonbondedMethod' : app.CutoffNonPeriodic})
```

To use the [Open Force Field `openff-1.2.0`](https://github.com/openforcefield/openff-forcefields), an update of the [Open Force Field ("Parsley") small molecule force field](https://openforcefield.org/news/introducing-openforcefield-1.0/) instead of GAFF 2.11, we would have instead specified `small_molecule_forcefield='openff-1.2.0'`.

To use [espaloma](https://github.com/choderalab/espaloma) for assigning small molecule parameters, for example with the [`espaloma-0.2.0` model](https://github.com/choderalab/espaloma/releases/tag/0.2.0) released with the [espaloma preprint](https://arxiv.org/abs/2010.01196), you can specify `small_molecule_forcefield='espaloma-0.2.0'`.

# Frequently Asked Questions (FAQ)

**Q:** What is the minimum version of OpenMM required to use this package?
<br>
**A:** You need at least OpenMM 7.4.2 to use the `openmmforcefields` package.

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
