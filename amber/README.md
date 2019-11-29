# Amber forcefield conversion utilities for OpenMM

## Manifest

* `convert_amber.py` - Amber forcefield conversion drive
* `biopolymer.yaml` - YAML file directing the conversion of AMBER biopolymer force fields
* `gaff.yaml` - YAML file directing conversion of GAFF small molecule force fields
* `merge_lipids.py` - Merge Amber split-residue lipids into single-residue OpenMM `ffxml` residue definitions
* `files/` - miscellaneous files used for testing and validation
* `ffxml/` - converted AMBER biopolymer force fields as OpenMM `ffxml`
* `gaff/dat` - historical versions of GAFF `.dat` force field parameter files
* `gaff/ffxml` - converted GAFF small molecule force fields as OpenMM `ffxml`
* `test/` - input files for testing conversion produces correct energies

## Amber forcefield conversion: `convert_amber.py`

Run `python convert_amber.py -h` to see all help options:
```
usage: convert_amber.py [-h] [--input INPUT] [--input-format INPUT_FORMAT]
                        [--verbose] [--log LOGFILE] [--protein-test] [--nucleic-test]
                        [--protein-ua-test] [--phospho-protein-test] [--gaff-test]

AMBER --> OpenMM forcefield conversion script

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        path of the input file. Default: "master.yaml"
  --input-format INPUT_FORMAT, -if INPUT_FORMAT
                        format of the input file: "yaml" or "leaprc". Default:
                        "yaml"
  --verbose, -v         turns verbosity on
  --log LOGFILE         log tests to specified CSV file
  --protein-test        validate resulting XML through protein tests
  --nucleic-test        validate resulting XML through nucleic acid tests
  --protein-ua-test     validate resulting XML through united-atom protein
                        tests
  --phospho-protein-test
                        validate resulting XML through phosphorylated protein
                        tests
  --gaff-test           validate resulting XML through small-molecule (GAFF)
                        test
```

### Converting the AMBER force fields

* Install the appropriate `AmberTools`
```bash
conda install --yes ambertools==19.9
```
* Convert biopolymer force fields
```bash
python convert_amber.py --input biopolymer.yaml --log biopolymer-tests.csv
```
* Convert GAFF small molecule force fields:
```bash
python convert_amber.py --input gaff.yaml --log gaff-tests.csv
```

With the defaults as set, all you need to do is have the script and the `files/` directory and call `python convert_amber.py`. `-v` or `--no-log` as wanted.

Output:
* `ffxml/` - converted AMBER biopolymer force fields as OpenMM `ffxml`
* `gaff/ffxml` - converted GAFF small molecule force fields as OpenMM `ffxml`
* `biopolymer-log.csv` - CSV log for biopolymer energy discrepancy tests
* `gaff-tests.csv` - CSV log for GAFF energy discrepancy tests
* (LeAP is called extensively by the script, and outputs its own `leap.log` too)

`-v` will turn on printing of the progress-tracking comments from the script. Warnings issued by ParmEd (warnings are by default filtered to error, but in some cases this had to be relaxed to 'always' because of the nature of the files converted) are printed always. The output of LeAP is always redirected to `dev/null`.

### YAML input format

By default the script takes a YAML file input, `files/master.yaml` has everything that is going on here.
There's only a few rules to the required structure of the YAML and it will be very easily extendable to future forcefields.

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
Water models are converted manually---we supply them as an ffxml in `files/`.
`tip3p.xml`, `tip4pew.xml`, and `spce.xml` are provided from the OpenMM 7.0, with some changes to make them the 'newest' format (e.g. no classes used, only types).
These `ffxml`s are integrated together with the converted ion parameters to make all the output.
OpenMM 7.0 is therefore listed as the source of these files - hence the extra input.
```yaml
- sourcePackage2: OpenMM
  sourcePackageVersion2: 7.0
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

* `Source`: AMBER input files

* `Solvent_source`: the water file in `files/` for the standard (i.e. water model containing) XMLs **or** Standard - this is the same as the `Name` field for the appropriate standard (water model containing) XML - we need to know that, because for the 'overloading' sets both that XML and the standard XML need to be loaded for energy testing

* `Solvent`: this is the name of the solvent, this is necessary to avoid hardcoding of recognition of what solvent you're using from the names of the files etc. - and knowing which solvent you're using is necessary for energy validations.

* `Name`: the desired name of the ffxml. (For proteins and nucleic this is done by the script, which a product of `leaprc.ff14SB` will call `ff14SB.xml` etc.)

Should you want to provide a different YAML (shorter version, different completely, whatever you want), script will take it with:
```bash
python convert_amber.py --input name_of_your_yaml.yaml
```
The outputs of any YAML will be written to `ffxml/`.

You can also provide a `leaprc` of your choosing via:
```bash
python convert_amber.py --input name_of_your_leaprc --input-format leaprc
```
The output of any leaprc will be written to the `./`.

**A few remarks on the water models and ions**

For this conversion:
* water models converted manually using `ffxml` files placed in `files/`, and merged with the appropriate converted ions
* we create *standard* recommended combinations of water and ion models: `tip3p_standard.xml`, `tip4pew_standard.xml`, `spce_standard.xml`: water model + JC monovalent ions + compromise set +2 ions
* for each water model, we have an HFE and IOD set for multivalent ions; all have templates set to `overload = "1"`. (`tip3p_HFE_multivalent.xml`, `tip3p_IOD_multivalent.xml` etc.)
* usage is to always load in a standard, and then you can overload +2's and add +3 and +4 with the HFE or IOD files
* naming of the water atom types remains as before (`tip3p-O`)
* naming of the ion atom types is `name_of_set (dash) amber_atom_type_name`, e.g. `tip3p_standard-Na+`, `tip3p_HFE_multivalent-Zn2+`.

**Tests**

Run `test/test.py` with `nosetests`.

## Merging Amber split-residue lipids into single residues for OpenMM: `merge_lipids.py`

```bash
python merge_lipids.py
```

**Acknowledgments**

* Rafal Wiewiora (MSKCC) for creating these tools
* Junmei Wang (University of Pittsburgh) for assistance in compiling historical GAFF releases
