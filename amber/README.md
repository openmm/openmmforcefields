# Amber forcefield conversion utilities for OpenMM

## Manifest

* `convert.sh` - Driver script to regenerate all force fields
* `convert_amber.py` - Main script for converting Amber force fields
* `convert_amber_ions.py` - Script for converting Amber ions
* `biopolymer.yaml` - YAML file directing the conversion of Amber biopolymer force fields
* `gaff.yaml` - YAML file directing conversion of GAFF small molecule force fields
* `glycam/glycan.yaml` - YAML file directing conversion of Amber GLYCAM force fields
* `solvents.yaml` - YAML file directing conversion of solvents
* `files/` - miscellaneous files used for testing and validation
* `test/` - input files for testing conversion produces correct energies

## Amber forcefield conversion: `convert_amber.py`

Run `python convert_amber.py -h` to see all help options:
```
usage: convert_amber.py [-h] [--input INPUT] [--input-format INPUT_FORMAT]
                        [--output-dir OUTPUT_DIR] [--verbose]
                        [--log LOG_FILENAME] [--protein-test] [--nucleic-test]
                        [--protein-ua-test] [--phospho-protein-test]
                        [--gaff-test] [--lipids-test] [--combination-tests]

AMBER --> OpenMM forcefield conversion script

options:
  -h, --help            show this help message and exit
  --input, -i INPUT     path of the input file. Default: "master.yaml"
  --input-format, -if INPUT_FORMAT
                        format of the input file: "yaml" or "leaprc". Default:
                        "yaml"
  --output-dir, -od OUTPUT_DIR
                        path of the output directory. Default: "ffxml/" for
                        yaml, "./" for leaprc
  --verbose, -v         turns verbosity on
  --log LOG_FILENAME    log energies for tests to specified CSV file
  --protein-test        validate resulting XML through protein tests
  --nucleic-test        validate resulting XML through nucleic acid tests
  --protein-ua-test     validate resulting XML through united-atom protein
                        tests
  --phospho-protein-test
                        validate resulting XML through phosphorylated protein
                        tests
  --gaff-test           validate resulting XML through small-molecule (GAFF)
                        test
  --lipids-test         validate resulting XML through lipids tests
  --combination-tests   validate combinations of force fields
```

## Converting the AMBER force fields

Install the appropriate `AmberTools` from `conda-forge`, and run the main
conversion script.  This calls `convert_amber.py` and `convert_amber_ions.py`.
```bash
conda install -c conda-forge --yes ambertools
./convert.sh
```
The outputs will be placed in `openmmforcefields/ffxml/amber`.

For more control, you can run `convert_amber.py` directly given a YAML file
specifying a conversion to perform (the syntax is described below):
```bash
python convert_amber.py --input name_of_your_yaml.yaml
```
You can also provide a `leaprc` of your choosing via:
```bash
python convert_amber.py --input name_of_your_leaprc --input-format leaprc
```

## YAML input format

By default the script takes a YAML file input.  There's only a few rules to the
required structure of the YAML and it will be very easily extendable to future
forcefields.

First entry in the YAML must be:
```yaml
- sourcePackage: AmberTools
  sourcePackageVersion: 15
```
a MODE declaration follows:
```yaml
- MODE: LEAPRC
```
There are two `MODE`s: `LEAPRC` and `RECIPE`.
* `LEAPRC`: Convert the contents of a `leaprc` file.
* `RECIPE`: Use a mix of `.dat`, `frcmod` and `.lib`, rather than `leaprc`; used for water-ion conversions

`LEAPRC` mode is used for all protein/nucleic acids force fields, e.g.:
```yaml
- Source: leaprc.ff14SB
  Reference:
  - >-
    Maier, J.A., Martinez, C., Kasavajhala, K., Wickstrom, L., Hauser, K.E., and Simmerling, C. (2015).
    ff14SB: Improving the Accuracy of Protein Side Chain and Backbone Parameters from ff99SB. J. Chem. Theory Comput. 11, 3696-3713.
  Test:
  - protein
  - nucleic
```
There's an optional `Options` field which allows changes to parameters if the default must be overridden:
```yaml
- Source: leaprc.phosaa10
  Reference:
  - >-
    Steinbrecher, T., Latzer, J., and Case, D.A. (2012). Revised AMBER parameters for bioorganic phosphates.
    J. Chem. Theory Comput. 8, 4405-4412.
  Options:
    filter_warnings: always
    write_unused: True
  Test:
  - protein_phospho
```
For converting water and ion force fields, the `MODE` is changed to `RECIPE`:
```yaml
- MODE: RECIPE
```
and an extra source package is declared.

### Notes on conversion for water and ions

Water models are converted manually: we supply them as an ffxml in `files/`.
`tip3p.xml`, `tip4pew.xml`, and `spce.xml` are provided from OpenMM 7.5.0, with
some changes to make them the 'newest' format (e.g. no classes used, only
types).  These FFXMLs are integrated together with the converted ion parameters
to make all the output.  OpenMM 7.5.0 is therefore listed as the source of these
files: hence the extra input.
```yaml
- sourcePackage2: OpenMM
  sourcePackageVersion2: 7.5.0
```

Here are two examples:
1. A 'standard' file - water model + JC monovalent ions + compromise set +2 ions
```yaml
- Source:
  - parm/frcmod.ionsjc_tip3p
  - parm/frcmod.ionslrcm_cm_tip3p
  - lib/atomic_ions.lib
  Solvent_source: tip3p.xml
  Solvent: tip3p
  Name: tip3p_standard
  Reference:
  - >-
    Joung, I.S., and Cheatham, Thomas E. (2008).
    Determination of Alkali and Halide Monovalent Ion Parameters for Use in Explicitly Solvated Biomolecular Simulations. J. Phys. Chem. B 112, 9020-9041.
  - >-
    Joung, I.S., and Cheatham, T.E. (2009).
    Molecular Dynamics Simulations of the Dynamic and Energetic Properties of Alkali and Halide Ions Using Water-Model-Specific Ion Parameters. J. Phys. Chem. B 113, 13279-13290.
  - >-
    Li, P., Roberts, B.P., Chakravorty, D.K., and Merz, K.M. (2013).
    Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent. J. Chem. Theory Comput. 9, 2733-2748.
  - >-
    Jorgensen, W.L., Chandrasekhar, J., Madura, J.D., Impey, R.W., and Klein, M.L. (1983).
    Comparison of simple potential functions for simulating liquid water. The Journal of Chemical Physics 79, 926-935.
  Test:
  - water_ion
```
2. An 'overloading' set: HFE +2, +3 and +4 ions for tip3p water.
```yaml
- Source:
  - parm/frcmod.ionslrcm_hfe_tip3p
  - parm/frcmod.ions34lsm_hfe_tip3p
  - lib/atomic_ions.lib
  Standard: tip3p_standard
  Solvent: tip3p
  Name: tip3p_HFE_multivalent
  Reference:
  - >-
    Li, P., Roberts, B.P., Chakravorty, D.K., and Merz, K.M. (2013).
    Rational Design of Particle Mesh Ewald Compatible Lennard-Jones Parameters for +2 Metal Cations in Explicit Solvent. J. Chem. Theory Comput. 9, 2733-2748.
  - >-
    Li, P., Song, L.F., and Merz, K.M. (2015).
    Parameterization of Highly Charged Metal Ions Using the 12-6-4 LJ-Type Nonbonded Model in Explicit Water. J. Phys. Chem. B 119, 883-895.
  Test:
  - water_ion
```

Notes on syntax:
* `Source`: AMBER input files
* `Solvent_source`: the water file in `files/` for the standard (i.e. water model containing) XMLs **or** Standard - this is the same as the `Name` field for the appropriate standard (water model containing) XML - we need to know that, because for the 'overloading' sets both that XML and the standard XML need to be loaded for energy testing
* `Solvent`: this is the name of the solvent, this is necessary to avoid hardcoding of recognition of what solvent you're using from the names of the files etc. - and knowing which solvent you're using is necessary for energy validations.
* `Name`: the desired name of the ffxml. (For proteins and nucleic this is done by the script, which a product of `leaprc.ff14SB` will call `ff14SB.xml` etc.)

Notes on conversion in general:
* water models converted manually using `ffxml` files placed in `files/`, and merged with the appropriate converted ions
* we create *standard* recommended combinations of water and ion models: `tip3p_standard.xml`, `tip4pew_standard.xml`, `spce_standard.xml`: water model + JC monovalent ions + compromise set +2 ions
* for each water model, we have an HFE and IOD set for multivalent ions; all have templates set to `overload = "1"`. (`tip3p_HFE_multivalent.xml`, `tip3p_IOD_multivalent.xml` etc.)
* usage is to always load in a standard, and then you can overload +2's and add +3 and +4 with the HFE or IOD files
* naming of the water atom types remains as before (`tip3p-O`)
* naming of the ion atom types is `name_of_set (dash) amber_atom_type_name`, e.g. `tip3p_standard-Na+`, `tip3p_HFE_multivalent-Zn2+`.

## Acknowledgments

* Rafal Wiewiora (MSKCC) for creating these tools
* Junmei Wang (University of Pittsburgh) for assistance in compiling historical GAFF releases
