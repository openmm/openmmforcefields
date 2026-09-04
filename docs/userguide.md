# User Guide

## Using the AMBER and CHARMM biopolymer force fields

This repository contains force fields for use with the [OpenMM `ForceField` class](http://docs.openmm.org/latest/userguide/application.html#force-fields) for parameterizing biomolecular systems.
If you're not familiar with this approach to applying parameters to biomolecular systems, please see the [OpenMM User Guide](http://docs.openmm.org/latest/userguide/application.html#force-fields).

### Using the AMBER force fields

Once installed, the AMBER force fields will be registered in the `amber/` relative path searched by [`openmm.app.ForceField`](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField).

For example, to specify the newer recommended [`ff14SB`](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255) force field and accompanying recommended ions and solvent models (corresponding to force fields loaded in LEaP with `leaprc.protein.ff14SB`), prepend the `amber` prefix and the `.xml` suffix:
```python
forcefield = ForceField("amber/protein.ff14SB.xml")
```
To access just the [`ff14SB`](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.5b00255) force field converted from `oldff/leaprc.ff14SB`:
```python
forcefield = ForceField("amber/ff14SB.xml")
```
or to specifically access the older (now outdated) `ff99SBildn` force field converted from `oldff/leaprc.ff14SB`:
```python
forcefield = ForceField("amber/ff99SBildn.xml")
```
The TIP3P conversion also includes the [Joung and Cheatham recommended salt models](https://doi.org/10.1021/jp8001614) (`parm/frcmod.ionsjc_tip3p`) and recommended divalent counterion parameters (`parm/frcmod.ions234lm_126_tip3p`):
```python
forcefield = ForceField(
    "amber/protein.ff14SB.xml",
    "amber/tip3p_standard.xml",
    "amber/tip3p_HFE_multivalent.xml",
)
```

### Using the CHARMM force fields

Similarly, the CHARMM force fields will be registered in the `charmm/` relative path:
```python
forcefield = ForceField("charmm/charmm36.xml")
```

## Using AMBER GAFF 1.x and 2.x for small molecules

The OpenMMForceFields package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using [GAFF](http://ambermd.org/antechamber/gaff.html) versions 1 or 2.

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

If the provided molecule(s) contain nonzero (user-specified) partial charges (stored in the `Molecule.partial_charges` attribute), those partial charges will be used.
If they do _not_ contain partial charges, the OpenFF Toolkit is used to generate them.

* In `GAFFTemplateGenerator`, [Antechamber](http://ambermd.org/antechamber/) from the [AmberTools](http://ambermd.org/AmberTools.php) distribution (which uses the `sqm` semiempirical quantum chemical package) is used to assign AM1-BCC charges (`antechamber -c bcc`). The conformers used in charge assignment are generated by RDKit.
* In `SMIRNOFFTemplateGenerator`, the charges are assigned by the OpenFF software using [OpenFF NAGL](https://github.com/openforcefield/openff-nagl) or according to the [AM1-BCC settings of the SMIRNOFF specification](https://openforcefield.github.io/standards/standards/smirnoff/#toolkitam1bcc-temporary-support-for-toolkit-based-am1-bcc-partial-charges) (depending on the force field in use).
  * If OpenEye Toolkits are installed and licensed, OpenEye's Quacpac Toolkit is used to generate [AM1-BCC-ELF10 charges](https://docs.eyesopen.com/toolkits/python/quacpactk/OEProtonClasses/OEAM1BCCELF10Charges.html). In this case, the conformers used to generate charges are generated by OpenEye's Omega Toolkit.
  * If OpenEye Toolkits are not installed or not licensed, AmberTools's Antechamber/`sqm` is used to generate AM1-BCC charges in the same way that it is used in `GAFFTemplateGenerator`. In this case, the conformers used to generate charges are generated by RDKit.

*Note:* The `Molecule` object must have the all protons and stereochemistry explicitly specified, and must match the exact protonation and tautomeric state of the molecule that will be found in your OpenMM `Topology` object.
The atom ordering need not be the same.

*Note:* The first time a `Molecule` is specified, added, or cached, if it lacks partial charges, the automatically generated charges will be cached and reused; if it contains user-specified partial charges, those charges will be used and cached. Adding the molecule again with a different set of charges will have no effect on changing which charges are assigned.

### Caching

`GAFFTemplateGenerator` also supports the ability to specify a cache filename, allowing parameters for small molecules to be computed only once and then cached in the specified cache file thereafter.

### Examples using `GAFFTemplateGenerator` to generate small molecule GAFF parameters

Create a GAFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit import Molecule

molecule = Molecule.from_smiles("c1ccccc1")
# Create the GAFF template generator
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
)

gaff = GAFFTemplateGenerator(molecules=molecule, forcefield="gaff-2.2.20")
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField

forcefield = ForceField(
    "amber/protein.ff14SB.xml",
    "amber/tip3p_standard.xml",
    "amber/tip3p_HFE_multivalent.xml",
)
# Register the GAFF template generator
forcefield.registerTemplateGenerator(gaff.generator)
# You can now parameterize an OpenMM Topology object that contains the specified molecule.
# forcefield will load the appropriate GAFF parameters when needed, and antechamber
# will be used to generate small molecule parameters on the fly.
from openmm.app import PDBFile

pdbfile = PDBFile("t4-lysozyme-L99A-with-benzene.pdb")
system = forcefield.createSystem(pdbfile.topology)
```

Create a template generator for multiple molecules read from an SDF file:

```python
molecules = Molecule.from_file("molecules.sdf")
gaff = GAFFTemplateGenerator(molecules=molecules, forcefield="gaff-2.2.20")
```
You can also add molecules to the generator later, even after the generator has been registered:
```python
gaff.add_molecules(molecule)
gaff.add_molecules([molecule1, molecule2])
```
To check which GAFF versions are supported, examine the `INSTALLED_FORCEFIELDS` attribute:
```pycon
>>> print(GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)
['gaff-1.4', 'gaff-1.8', 'gaff-1.81', 'gaff-2.1', 'gaff-2.11', 'gaff-2.2.20']
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
generator = GAFFTemplateGenerator(cache="gaff-molecules.json", forcefield="gaff-1.8")
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Using the Open Force Field Initiative SMIRNOFF small molecule force fields

The OpenMMForceFields package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that automatically generates OpenMM residue templates for small molecules lacking parameters using the [Open Force Field Initiative](http://openforcefield.org) [SMIRNOFF](https://openforcefield.github.io/standards/standards/smirnoff/) small molecule force fields.
This includes the [`openff-1.x.y` ("Parsley")](https://openforcefield.org/news/introducing-openforcefield-1.0/) and [`openff-2.x.y` ("Sage")](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00039) small molecule force field lines, including the [most recent force field in each lines](https://github.com/openforcefield/openff-forcefields).

The `SMIRNOFFTemplateGenerator` residue template generator operates in a manner very similar to `GAFFTemplateGenerator`, so we only highlight its differences here.

### Examples using `SMIRNOFFTemplateGenerator` to generate small molecule SMIRNOFF parameters

Create a SMIRNOFF template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:

```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit import Molecule

molecule = Molecule.from_smiles("c1ccccc1")
# Create the SMIRNOFF template generator with the openff-2.3.0 force field
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield="openff-2.3.0")
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField

forcefield = ForceField(
    "amber/protein.ff14SB.xml",
    "amber/tip3p_standard.xml",
    "amber/tip3p_HFE_multivalent.xml",
)
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(smirnoff.generator)

# create a System with the non-bonded settings of mainline OpenFF force fields
# (9 Angstrom cut-off, switching distance applied at 8 Angstrom)
import openmm.unit
system = forcefield.createSystem(
    topology=molecule.to_topology().to_openmm(),
    nonbondedCutoff=0.9 * openmm.unit.nanometer,
    switchDistance=0.8 * openmm.unit.nanometer,
)
```

You can create a template generator capable of recognizing multiple molecules,
*e.g.*, read from an SDF file:

```python
molecules = Molecule.from_file("molecules.sdf")
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield="openff-2.3.0")
```

You can also add molecules to the generator later, even after the generator has been registered:

```python
smirnoff.add_molecules(molecule)
smirnoff.add_molecules([molecule1, molecule2])
```

You must specify the desired SMIRNOFF force field with the `forcefield` keyword
argument when you create a `SMIRNOFFTemplateGenerator` by providing at least one
force field name, path to an OFFXML force field file, file object, or XML in a
string.  An iterable can be provided to the `forcefield` argument to load from
multiple sources if desired.  For exapmle:

```python
# Create a SMIRNOFF residue template generator from an older Sage or Parsley version
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield="openff-1.3.0")
# Request a specific SMIRNOFF force field installed with OpenFF by its full name
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield="openff_unconstrained-1.3.0.offxml")
# Use a local .offxml file instead
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecules, forcefield="local-file.offxml")
# Load from multiple sources at once
smirnoff = SMIRNOFFTemplateGenerator(molecules=molecule, forcefield=['openff-2.3.0', 'tip5p.offxml'])
```

You can check the full paths of the force field files that have been loaded:

```pycon
>>> smirnoff.smirnoff_filenames
['/.../offxml/openff_unconstrained-2.3.0.offxml', '/.../offxml/tip5p.offxml']
```

Any force field installed with the OpenFF Toolkit can be requested by its
filename, *e.g.*, `openff_unconstrained-2.3.0.offxml`.  To see which force
fields can be requested by an abbreviated name, examine `INSTALLED_FORCEFIELDS`:

```pycon
>>> SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS
['openff-1.0.0', 'openff-1.0.1', 'openff-1.1.0', 'openff-1.1.1', ...]
```

The names in this list are of the form `openff-x.y.z` and will resolve to the
"unconstrained" variants of standard OpenFF force fields, which can also be
loaded with filenames of the form `openff_unconstrained-x.y.z.offxml`.  To
obtain the "constrained" variants, use filenames `openff-x.y.z.offxml`.  When
the constrained variants are used, bonds to hydrogen atoms in molecules
parameterized by the template generator will be constrained unconditionally,
regardless of the constraint options given in `ForceField.createSystem()`.  The
unconstrained variants, which are the defaults selected if one of the predefined
names in `INSTALLED_FORCEFIELDS` is used, do not force these constraints to be
present, and the option of which degrees of freedom should be constrained are
left up to the user by passing any desired `constraint` setting to
[`ForceField.createSystem()`](https://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField.createSystem).

You can optionally specify a file that contains a cache of pre-parameterized molecules:

```python
smirnoff = SMIRNOFFTemplateGenerator(
    cache="smirnoff-molecules.json",
    forcefield="openff-2.3.0",
)
```

Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Using espaloma to generate small molecule force fields

The OpenMMForceFields package includes a [residue template generator](http://docs.openmm.org/latest/userguide/application.html#adding-residue-template-generators) for [the OpenMM `ForceField` class](http://docs.openmm.org/latest/api-python/generated/openmm.app.forcefield.ForceField.html#openmm.app.forcefield.ForceField) that can automatically generate OpenMM residue templates for small molecules lacking parameters using [espaloma](https://github.com/choderalab/espaloma) via one of its released force fields, provided `espaloma` and its dependencies are installed.
`espaloma` uses a [graph convolutional model](https://arxiv.org/abs/2010.01196) to generate both valence parameters and fast partial charges.


The `EspalomaTemplateGenerator` residue template generator operates in a manner very similar to `GAFFTemplateGenerator`, so we only highlight its differences here.

**This feature is currently experimental and its API is subject to change.**

### Examples using `EspalomaTemplateGenerator` to generate small molecule parameters with espaloma

Create an espaloma template generator for a single molecule (benzene, created from SMILES) and register it with ForceField:
```python
# Create an OpenFF Molecule object for benzene from SMILES
from openff.toolkit import Molecule

molecule = Molecule.from_smiles("c1ccccc1")
# Create the SMIRNOFF template generator with the released espaloma-0.3.2 force field
from openmmforcefields.generators import (
    EspalomaTemplateGenerator,
)

espaloma = EspalomaTemplateGenerator(molecules=molecule, forcefield="espaloma-0.3.2")
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from openmm.app import ForceField

forcefield = ForceField(
    "amber/protein.ff14SB.xml",
    "amber/tip3p_standard.xml",
    "amber/tip3p_HFE_multivalent.xml",
)
# Register the SMIRNOFF template generator
forcefield.registerTemplateGenerator(espaloma.generator)
```
Create a template generator for a specific espaloma force field for multiple molecules read from an SDF file:
```python
molecules = Molecule.from_file("molecules.sdf")
# Create an espaloma residue template generator from the espaloma-0.3.2 release,
# retrieving it automatically from GitHub release artifacts
espaloma = EspalomaTemplateGenerator(molecules=molecules, forcefield="espaloma-0.3.2")
# Create an espaloma residue template generator from an espaloma model retrieved from a URL
espaloma = EspalomaTemplateGenerator(
    molecules=molecules,
    forcefield="https://github.com/choderalab/espaloma/releases/download/0.3.2/espaloma-0.3.2.pt",
)
# Create an espaloma residue template generator from an espaloma model stored in a local file
espaloma = EspalomaTemplateGenerator(
    molecules=molecules,
    forcefield="/path/to/espaloma-0.3.2.pt",
)
```
You can also add molecules to the generator later, even after the generator has been registered:
```python
smirnoff.add_molecules(molecule)
smirnoff.add_molecules([molecule1, molecule2])
```
You can optionally specify a file that contains a cache of pre-parameterized molecules:
```python
espaloma = EspalomaTemplateGenerator(
    cache="espaloma-molecules.json",
    forcefield="espaloma-0.3.2",
)
```
Newly parameterized molecules will be written to the cache, saving time next time these molecules are encountered.

## Automating force field management with `SystemGenerator`

The OpenMMForceFields package provides the `openmmforcefields.generators.SystemGenerator` class that handles management of common force fields transparently for you.

### Using `SystemGenerator` to automate the use of AMBER force fields with GAFF, OpenFF, or espaloma for small molecule parameterization

Here's an example that uses GAFF 2.2.20 along with the new `ff14SB` generation of AMBER force fields (and compatible solvent models) to generate an OpenMM `System` object from an [Open Force Field `Topology`](https://open-forcefield-toolkit.readthedocs.io/en/latest/api/generated/openff.toolkit.topology.Topology.html#openff.toolkit.topology.Topology) object:
```python
# Define the keyword arguments to feed to ForceField
from openmm import unit
from openmm import app

forcefield_kwargs = {
    "constraints": app.HBonds,
    "rigidWater": True,
    "removeCMMotion": False,
    "hydrogenMass": 4 * unit.amu,
}
# Initialize a SystemGenerator using GAFF
from openmmforcefields.generators import SystemGenerator

system_generator = SystemGenerator(
    forcefields=[
        "amber/ff14SB.xml",
        "amber/tip3p_standard.xml",
    ],
    small_molecule_forcefield="gaff-2.2.20",
    forcefield_kwargs=forcefield_kwargs,
    cache="db.json",
)
# Create an OpenMM System from an OpenMM Topology object
system = system_generator.create_system(openmm_topology)
# Alternatively, create an OpenMM System from an OpenMM Topology object and a list of OpenFF Molecule objects
molecules = Molecule.from_file("molecules.sdf", file_format="sdf")
system = system_generator.create_system(openmm_topology, molecules=molecules)
```
Parameterized molecules are cached in `db.json`.
Parameters for multiple force fields can be held in the same cache file.

By default, `SystemGenerator` will use `PME` for periodic systems and `NoCutoff` for non-periodic systems.
You can modify this behavior with the optional `periodic_forcefield_kwargs` and `nonperiodic_forcefield_kwargs` arguments, which are used to update `forcefield_kwargs` depending on whether the system is periodic or non-periodic:
```python
from openmm import app

system_generator = SystemGenerator(
    forcefields=[
        "amber/ff14SB.xml",
        "amber/tip3p_standard.xml",
    ],
    periodic_forcefield_kwargs={"nonbondedMethod": app.LJPME},
    nonperiodic_forcefield_kwargs={"nonbondedMethod": app.CutoffNonPeriodic},
)
```

To use the [OpenFF's Sage `openff-2.2.1`](https://github.com/openforcefield/openff-forcefields) or a newer version, an update of the [Open Force Field ("Parsley") small molecule force field](https://openforcefield.org/news/introducing-openforcefield-1.0/) instead of GAFF 2.2.20, we would have instead specified `small_molecule_forcefield='openff-2.2.1'`.

To use [espaloma](https://github.com/choderalab/espaloma) for assigning small molecule parameters, for example with the [`espaloma-0.3.2` model](https://github.com/choderalab/espaloma/releases/tag/0.3.2) released with the [espaloma preprint](https://arxiv.org/abs/2307.07085), you can specify `small_molecule_forcefield='espaloma-0.3.2'`.
