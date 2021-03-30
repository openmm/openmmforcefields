# AMBER --> OpenMM force-field conversion script
# Author: Rafal P. Wiewiora, ChoderaLab
from __future__ import print_function, division
from io import StringIO
import parmed
from simtk import openmm
import simtk.openmm.app as app
import simtk.unit as u
import simtk
import os
import sys
import re
import tempfile
import yaml
from distutils.spawn import find_executable
import hashlib
from collections import OrderedDict
import glob
import argparse
from lxml import etree as et
import csv
import logging
import warnings
import xml.etree.ElementTree as etree
from copy import deepcopy
from parmed.exceptions import ParameterWarning
warnings.filterwarnings('error', category=ParameterWarning)

_loadoffre = re.compile(r'loadoff (\S*)', re.I)
_sourcere = re.compile(r'source (\S*)', re.I)

# check for AMBERHOME, find from tleap location if not set, exception if can't
if os.getenv('AMBERHOME'):
    AMBERHOME = os.getenv('AMBERHOME')
else:
    if not find_executable('tleap'):
        raise Exception('AMBERHOME not set and tleap not available from PATH')
    tleap_path = find_executable('tleap')
    AMBERHOME = os.path.split(tleap_path)[0]
    AMBERHOME = os.path.join(AMBERHOME, '../')
    parmed.amber.AMBERHOME = AMBERHOME

# set global defaults for verbose and log
verbose = False
no_log = False

# set files that are ignored in leaprc's
# solvents and ions converted separately; leaprc.ff10 calls phosphoaa10.lib
# which does not exist anymore, LeAP skips it on error so we do too
ignore = {'solvents.lib', 'atomic_ions.lib', 'ions94.lib', 'ions91.lib',
          'phosphoaa10.lib'}

# define NEARLYZERO to replace numerical comparisons to zero
NEARLYZERO = 1e-10

# set beta
temperature = 300.0 * u.kelvin
kB = u.BOLTZMANN_CONSTANT_kB * u.AVOGADRO_CONSTANT_NA
kT = kB * temperature
beta = 1.0/kT

class LeapException(Exception):
    def __init__(self, leaprc_filename):
        msg = 'Something went wrong in processing this LEaP input file:\n'
        msg += '\n'
        infile = open(leaprc_filename, 'rt')
        contents = infile.read()
        msg += contents
        msg += '\n'
        super(LeapException, self).__init__(msg)

def main():
    global verbose
    global no_log
    global logger
    # argparse
    parser = argparse.ArgumentParser(description='AMBER --> OpenMM forcefield '
                                                 'conversion script')
    parser.add_argument('--input', '-i', default='master.yaml',
                        help='path of the input file. Default: "master.yaml"')
    parser.add_argument('--input-format', '-if', default='yaml',
                        help='format of the input file: "yaml" or "leaprc". Default: "yaml"')
    parser.add_argument('--output-dir', '-od', help='path of the output directory. '
                        'Default: "ffxml/" for yaml, "./" for leaprc')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='turns verbosity on')
    parser.add_argument('--log', action='store', dest='log_filename', default=None,
                        help='log energies for tests to specified CSV file')
    parser.add_argument('--protein-test', action='store_true',
                        help='validate resulting XML through protein tests')
    parser.add_argument('--nucleic-test', action='store_true',
                        help='validate resulting XML through nucleic acid tests')
    parser.add_argument('--protein-ua-test', action='store_true',
                        help='validate resulting XML through united-atom protein tests')
    parser.add_argument('--phospho-protein-test', action='store_true',
                        help='validate resulting XML through phosphorylated protein tests')
    parser.add_argument('--gaff-test', action='store_true',
                        help='validate resulting XML through small-molecule (GAFF) test')
    parser.add_argument('--lipids-test', action='store_true',
                        help='validate resulting XML through lipids tests')
    args = parser.parse_args()
    verbose = args.verbose

    if args.log_filename:
        logger = Logger(args.log_filename) # log to file
    else:
        logger = Logger() # be silent

    # input is either a YAML or a leaprc - default is leaprc
    # output directory hardcoded here for ffxml/
    if args.input_format == 'yaml':
        if args.output_dir is None:
            convert_yaml(args.input, ffxml_dir='ffxml/')
        else:
            convert_yaml(args.input, ffxml_dir=args.output_dir)
    # if leaprc converted - output to the same dir
    elif args.input_format == 'leaprc':
        if args.output_dir is None:
            ffxml_name = convert_leaprc(args.input, ffxml_dir='./')
        else:
            ffxml_name = convert_leaprc(args.input, ffxml_dir=args.output_dir)
        if args.protein_test:
            validate_protein(ffxml_name, args.input)
        if args.nucleic_test:
            validate_nucleic(ffxml_name, args.input)
        if args.protein_ua_test:
            validate_protein(ffxml_name, args.input, united_atom=True)
        if args.phospho_protein_test:
            validate_phospho_protein(ffxml_name, args.input)
        if args.gaff_test:
            validate_gaff(ffxml_name, args.input)
        if args.lipids_test:
            validate_lipids(ffxml_name, args.input)
    else:
        sys.exit('Wrong input_format chosen.')

    logger.close()

def read_lines(filename):
    """
    Read lines from file, stripping comments and newlines.

    """
    with open(filename, 'rt') as f:
        lines = [ line if '#' not in line else line[:line.index('#')] for line in f.readlines() ]
    return lines

def write_file(file, contents):
    """
    Write text to file.

    Parameters
    ----------
    filename : str
       The file to write to
    contents : str
       Text contents to be written to file

    """
    if type(file) == str:
        outfile = open(file, 'w')
    else:
        outfile = os.fdopen(file, 'w')
    outfile.write(contents)
    outfile.close()

def convert_leaprc(files, split_filename=False, ffxml_dir='./', ignore=ignore,
    provenance=None, write_unused=False, override_level=None, filter_warnings='error'):
    if verbose: print('\nConverting %s to ffxml...' % files)
    # allow for multiple source files - further code assuming list is passed
    if not isinstance(files, list):
        files = [files]
    basename = ''
    for f in files:
        f_basename = os.path.basename(f)
        if split_filename:
            f_basename = f_basename.split('.')[1:]
            f_basename = '.'.join(f_basename)
        if not basename:
            basename = f_basename
        else:
            basename += '_'
            basename += f_basename
    ffxml_name = os.path.join(ffxml_dir, (basename + '.xml'))
    if not os.path.exists(ffxml_dir):
        os.mkdir(ffxml_dir)
    if verbose: print('Preprocessing the leaprc for %s...' % basename)
    # do source processing
    new_files = []
    for fil in files:
        lines = read_lines(fil)
        for line in lines:
            if _sourcere.findall(line):
                replace_leaprc = _sourcere.findall(line)[0]
                replace_leaprc_path = os.path.join(os.path.join(AMBERHOME,
                'dat/leap/cmd', replace_leaprc))
                new_files.append(replace_leaprc_path)
        new_files.append(fil)
    # now do ignore processing and join multiple leaprc's
    files = new_files
    new_lines = []
    for fil in files:
        lines = read_lines(fil)
        fil_new_lines = []
        for line in lines:
            if (ignore is not None and _loadoffre.findall(line) and
            _loadoffre.findall(line)[0] in ignore):
                continue
            fil_new_lines += line
        new_lines += fil_new_lines
    leaprc = StringIO(''.join(new_lines))
    if verbose: print('Converting to ffxml %s...' % ffxml_name)
    print("leaprc", ''.join(new_lines))
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, remediate_residues=(not write_unused))
    if override_level:
        for name, residue in params.residues.items():
            residue.override_level = override_level
    if filter_warnings != 'error':
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_name, provenance=provenance, write_unused=write_unused, improper_dihedrals_ordering='amber')
    else:
        params.write(ffxml_name, provenance=provenance, write_unused=write_unused, improper_dihedrals_ordering='amber')
    if verbose: print('%s successfully written!' % ffxml_name)
    return ffxml_name

def convert_gaff(files, ffxml_basename='', split_filename=False, ffxml_dir='./', ignore=ignore,
    provenance=None, write_unused=False, filter_warnings='error'):
    if verbose: print('\nConverting %s to ffxml...' % files)
    # allow for multiple source files - further code assuming list is passed
    if not isinstance(files, list):
        files = [files]
    # Create ffxml
    ffxml_name = os.path.join(ffxml_dir, (ffxml_basename + '.xml'))
    if not os.path.exists(ffxml_dir):
        os.mkdir(ffxml_dir)
    # Process parameter file
    params = parmed.amber.AmberParameterSet(files)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params, remediate_residues=(not write_unused))
    if filter_warnings != 'error':
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_name, provenance=provenance, write_unused=write_unused, improper_dihedrals_ordering='amber')
    else:
        params.write(ffxml_name, provenance=provenance, write_unused=write_unused, improper_dihedrals_ordering='amber')
    if verbose: print('%s successfully written!' % ffxml_name)
    return ffxml_name

def convert_recipe(files, solvent_file=None, ffxml_dir='./', provenance=None, ffxml_basename=None,
                   filter_warnings='always'):
    if verbose: print('\nConverting %s to ffxml...' % files)
    ffxml_name = os.path.join(ffxml_dir, (ffxml_basename + '.xml'))
    ffxml_temp_stringio = StringIO()
    params = parmed.amber.AmberParameterSet(files)
    print(params.atom_types.keys())
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    # Change atom type naming
    # atom_types
    new_atom_types = OrderedDict()
    for name, atom_type in params.atom_types.items():
        new_name = ffxml_basename + '-' + name
        new_atom_types[new_name] = atom_type
    params.atom_types = new_atom_types
    # atoms in residues
    for name, residue in params.residues.items():
        for atom in residue:
            new_type = ffxml_basename + '-' + atom.type
            atom.type = new_type
    if solvent_file is None:
    # this means this file does not include a water model - hard-coded assumption it is
    # then a 'multivalent' file - set overrideLevel to 1 for all residue templates
        for name, residue in params.residues.items():
            residue.override_level = 1
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_name, provenance=provenance, write_unused=False, improper_dihedrals_ordering='amber')
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_temp_stringio, provenance=provenance, write_unused=False, improper_dihedrals_ordering='amber')
        ffxml_temp_stringio.seek(0)
        if verbose: print('Modifying converted ffxml to append solvent parameters')
        tree_main = et.parse(ffxml_temp_stringio)
        tree_water = et.parse(solvent_file)
        root_main = tree_main.getroot()
        root_water = tree_water.getroot()
        with open(ffxml_name, 'wb') as f:
            f.write(b'<ForceField>\n ')
            f.write(et.tostring(root_main.findall('Info')[0]))
            f.write(b'<AtomTypes>\n  ')
            for subelement in root_main.findall('AtomTypes')[0]:
                f.write(et.tostring(subelement))
            f.write(b' ')
            for subelement in root_water.findall('AtomTypes')[0]:
                f.write(et.tostring(subelement))
            f.write(b'</AtomTypes>\n <Residues>\n  ')
            for subelement in root_main.findall('Residues')[0]:
                f.write(et.tostring(subelement))
            f.write(b' ')
            for subelement in root_water.findall('Residues')[0]:
                f.write(et.tostring(subelement))
            f.write(b'</Residues>\n <HarmonicBondForce>\n  ')
            for subelement in root_water.findall('HarmonicBondForce')[0]:
                f.write(et.tostring(subelement))
            f.write(b'</HarmonicBondForce>\n <HarmonicAngleForce>\n  ')
            for subelement in root_water.findall('HarmonicAngleForce')[0]:
                f.write(et.tostring(subelement))
            f.write(b'</HarmonicAngleForce>\n ')
            f.write(('<NonbondedForce coulomb14scale="%s" lj14scale="%s">\n  ' %
                   (root_main.findall('NonbondedForce')[0].attrib['coulomb14scale'],
                    root_main.findall('NonbondedForce')[0].attrib['lj14scale'])
                   ).encode('utf-8'))
            for subelement in root_main.findall('NonbondedForce')[0]:
                f.write(et.tostring(subelement))
            f.write(b' ')
            for subelement in root_water.findall('NonbondedForce')[0]:
                if subelement.tag == 'UseAttributeFromResidue': continue
                f.write(et.tostring(subelement))
            f.write(b'</NonbondedForce>\n</ForceField>')
    if verbose: print('%s successfully written!' % ffxml_name)
    return ffxml_name

def convert_yaml(yaml_name, ffxml_dir, ignore=ignore):
    data = yaml.load(open(yaml_name), Loader=yaml.FullLoader)
    # TODO: Verify that the version that is installed via conda matches sourcePackageVersion

    # Default yaml reading mode is leaprc
    ALLOWED_MODES = ('LEAPRC', 'RECIPE', 'GAFF')
    for entry in data:
        # Handle MODE switching
        if 'MODE' in entry:
            MODE = entry['MODE']
            if not MODE in ALLOWED_MODES:
                raise Exception(f'MODE definition must be one of {ALLOWED_MODES}')
            continue

        # Handle definition of source packages
        if 'sourcePackage' in entry:
            source_pack = entry['sourcePackage']
            source_pack_ver = entry['sourcePackageVersion']
            continue
        if 'sourcePackage2' in entry:
            # Switch mode to RECIPE processing
            source_pack2 = entry['sourcePackage2']
            source_pack_ver2 = entry['sourcePackageVersion2']
            continue

        # Extract source files, reference, and test files
        source_files = entry['Source']
        reference = entry['Reference']
        test_filename = entry['Test']

        # Make sure source_files is a list
        if isinstance(source_files, str):
            source_files = [source_files]

        # Recipes require extra definitions
        if MODE == 'RECIPE':
            recipe_name = entry['Name']
            solvent_name = entry['Solvent']
            if 'Solvent_source' in entry:
                recipe_source2 = entry['Solvent_source']
            else:
                recipe_source2 = None
            if 'Standard' in entry:
                standard_ffxml = os.path.join(ffxml_dir, (entry['Standard'] + '.xml'))
            else:
                standard_ffxml = None
        elif MODE == 'GAFF':
            recipe_name = entry['Name']

        # Create provenance object
        provenance = OrderedDict()
        files = []
        source = provenance['Source'] = []
        for source_file in source_files:
            if MODE == 'LEAPRC':
#                _filename = os.path.join(AMBERHOME, 'dat/leap/cmd', source_file)
                 _filename = os.path.join('./', source_file)
            elif MODE == 'RECIPE':
                _filename = os.path.join(AMBERHOME, 'dat/leap/', source_file)
            elif MODE == 'GAFF':
                _filename = os.path.join('gaff', 'dat', source_file)
            files.append(_filename)
            source.append(OrderedDict())
            source[-1]['Source'] = source_file
            md5 = hashlib.md5()
            with open(_filename, 'rb') as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]['md5hash'] = md5
            source[-1]['sourcePackage'] = source_pack
            source[-1]['sourcePackageVersion'] = source_pack_ver

        # For recipes, add water file and source info for it
        if MODE == 'RECIPE' and recipe_source2 is not None:
            _filename = os.path.join('files', recipe_source2)
            solvent_file = _filename
            source.append(OrderedDict())
            source[-1]['Source'] = recipe_source2
            md5 = hashlib.md5()
            with open(_filename, 'rb') as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]['md5hash'] = md5
            source[-1]['sourcePackage'] = source_pack2
            source[-1]['sourcePackageVersion'] = source_pack_ver2
        elif MODE == 'RECIPE' and recipe_source2 is None:
            solvent_file = None
        provenance['Reference'] = reference

        # set default conversion options
        write_unused = False
        filter_warnings = 'error'
        override_level = None
        # set conversion options if present
        if 'Options' in entry:
            for option in entry['Options']:
                if option == 'write_unused':
                    write_unused = entry['Options'][option]
                elif option == 'filter_warnings':
                    filter_warnings = entry['Options'][option]
                elif option == 'ffxml_dir':
                    ffxml_dir = entry['Options'][option]
                elif option == 'override_level':
                    override_level = entry['Options'][option]
                else:
                    raise Exception("Wrong option used in Options for %s"
                                       % source_files)

        # Convert files
        if MODE == 'LEAPRC':
            ffxml_name = convert_leaprc(files, ffxml_dir=ffxml_dir, ignore=ignore,
                         provenance=provenance, write_unused=write_unused, override_level=override_level,
                         filter_warnings=filter_warnings, split_filename=True)
        elif MODE == 'RECIPE':
            ffxml_name = convert_recipe(files, solvent_file=solvent_file,
                         ffxml_dir=ffxml_dir, provenance=provenance,
                         ffxml_basename=recipe_name)
        elif MODE == 'GAFF':
            ffxml_name = convert_gaff(files, ffxml_basename=recipe_name, ffxml_dir=ffxml_dir, ignore=ignore,
                         provenance=provenance, write_unused=write_unused,
                         filter_warnings=filter_warnings, split_filename=True)

        if 'CharmmFFXMLFilename' in entry:
            charmm_ffxml_filename = entry['CharmmFFXMLFilename']
            charmm_lipid2amber_filename = entry['CharmmLipid2AmberFilename']
            if verbose: print('Merging lipid entries...')
            merge_lipids(ffxml_name, charmm_ffxml_filename, charmm_lipid2amber_filename)
        if 'Prefix' in entry:
            prefix = entry['Prefix']
            if verbose: print('Rewriting %s to append prefix "%s"...' % (ffxml_name, prefix))
            add_prefix_to_ffxml(ffxml_name, prefix)
        if verbose: print('Validating the conversion...')
        tested = False
        for test in test_filename:
            if test == 'protein':
                validate_protein(ffxml_name, entry['Source'])
                tested = True
            elif test == 'nucleic':
                validate_dna(ffxml_name, entry['Source'])
                validate_rna(ffxml_name, entry['Source'])
                tested = True
            elif test == 'protein_ua':
                validate_protein(ffxml_name, entry['Source'], united_atom=True)
                tested = True
            elif test == 'protein_phospho':
                validate_phospho_protein(ffxml_name, entry['Source'])
                tested = True
            elif test == 'gaff':
                validate_gaff(ffxml_name, entry['leaprc'], entry['Source'])
                tested = True
            elif test == 'water_ion':
                validate_water_ion(ffxml_name, files, solvent_name, recipe_name,
                standard_ffxml=standard_ffxml)
                tested = True
            elif test == 'dna':
                validate_dna(ffxml_name, entry['Source'])
                tested = True
            elif test == 'rna':
                validate_rna(ffxml_name, entry['Source'])
                tested = True
            elif test == 'lipids':
                #validate_lipids(ffxml_name, source_files)
                validate_merged_lipids(ffxml_name, entry['Source'])
                tested = True
            elif test == 'protein_glycan':
                validate_glyco_protein(ffxml_name, entry['Source'])
                tested = True
        if not tested:
            raise Exception('No validation tests have been run for %s' %
                                source_files)

def merge_lipids(ffxml_filename, charmm_ffxml_filename, charmm_lipid2amber_filename):
    """
    Merge lipid residue definitions in AMBER ffxml file according to entries in a CHARMM ffxml file.

    Parameters
    ----------
    ffxml_filename : str
       AMBER lipids ffxml filename with AMBER split lipids.
    charmm_ffxml_filename : str
       CHARMM ffxml lipids
    charmmlipid2amber_filename : str
       CHARMM CSV file containing translation from CHARMM -> AMBER

    """
    # Read the input files.
    charmmff = etree.parse(charmm_ffxml_filename)
    amberff = etree.parse(ffxml_filename)
    charmmResidues = charmmff.getroot().find('Residues').findall('Residue')
    amberResidues = amberff.getroot().find('Residues').findall('Residue')
    amberResMap = {}
    for res in amberResidues:
        atoms = dict((atom.attrib['name'], atom) for atom in res.findall('Atom'))
        amberResMap[res.attrib['name']] = atoms
    translations = {}
    with open(charmm_lipid2amber_filename) as input:
        # Skip the first two lines.
        input.readline()
        input.readline()
        for line in input:
            fields = line.split(',')
            mergedRes = fields[0]
            mergedAtom = fields[2].split()[0]
            originalAtom, originalRes = fields[3].split()
            translations[(mergedRes, mergedAtom)] = (originalRes, originalAtom)

    # Remove all residues from the Amber file.
    parentNode = amberff.getroot().find('Residues')
    for res in amberResidues:
        parentNode.remove(res)

    # Copy over the CHARMM residues, making appropriate replacements.
    def translateResidue(residue):
        newres = deepcopy(residue)
        # Translate atom properties
        for atom in newres.findall('Atom'):
            key = (residue.attrib['name'], atom.attrib['name'])
            if key not in translations:
                return None # We don't have a translation.
            amberResName, amberAtomName = translations[key]
            if amberResName not in amberResMap or amberAtomName not in amberResMap[amberResName]:
                return None # We don't have a translation.
            amberAtom = amberResMap[amberResName][amberAtomName]
            for attrib in amberAtom.attrib:
                if attrib != 'name':
                    atom.attrib[attrib] = amberAtom.attrib[attrib]
        # Remove Patches from CHARMM residues
        for patch in newres.findall('AllowPatch'):
            newres.remove(patch)

        return newres

    # Iterate over CHARMM lipid residue templates and replace components with AMBER parameters
    for residue in charmmResidues:
        copy = translateResidue(residue)
        if copy is not None:
            parentNode.append(copy)

    # Write merged lipid ffxml file (overwriting original file)
    amberff.write(ffxml_filename)

def add_prefix_to_ffxml(ffxml_filename, prefix):
    """
    Replace the contents of an ffxml file with a modified version in which every atom type is prefixed with `prefix`.

    Parameters
    ----------
    ffxml_filename : str
       OpenMM ffxml filename (will be overwritten)
    prefix : str
       Prefix

    """

    import re
    import sys

    inTypes = False
    replacements = {}

    modified_contents = ''
    with open(ffxml_filename, 'r') as infile:
        for line in infile:
            if '<AtomTypes>' in line:
                inTypes = True
            if '</AtomTypes>' in line:
                inTypes = False
            if inTypes:
                match = re.search('name="(.*?)"', line)
                if match is not None:
                    name = match.group(1)
                    newName = prefix + '-' + name
                    line = line.replace('name="%s"' % name, 'name="%s"' % newName)
                    replacements['type="%s"' % name] = 'type="%s"' % newName
                    replacements['type1="%s"' % name] = 'type1="%s"' % newName
                    replacements['type2="%s"' % name] = 'type2="%s"' % newName
                    replacements['type3="%s"' % name] = 'type3="%s"' % newName
                    replacements['type4="%s"' % name] = 'type4="%s"' % newName
            else:
                for key in replacements:
                    if key in line:
                        line = line.replace(key, replacements[key])
            if line.endswith('\n'):
                line = line[:-1]
            modified_contents += line + '\n'

    with open(ffxml_filename, 'w') as outfile:
        outfile.write(modified_contents)

def assert_energies2(prmtop, inpcrd, ffxml, system_name='unknown', tolerance=2.5e-5,
                    improper_tolerance=1e-2, units=u.kilojoules_per_mole, openmm_topology=None, openmm_positions=None):

    # Get AMBER system
    parm_amber = parmed.load_file(prmtop, inpcrd)
    system_amber = parm_amber.createSystem()

    # Get OpenMM system
    if isinstance(ffxml, str):
        ff = app.ForceField(ffxml)
    else:
        ff = app.ForceField(*ffxml)

    if openmm_positions is None:
        openmm_positions = parm_amber.positions

    if openmm_topology is not None:
        system_omm = ff.createSystem(openmm_topology)
    else:
        system_omm = ff.createSystem(parm_amber.topology)

    def compute_potential_components(system, positions, beta=beta):
        # Note: this is copied from perses
        # TODO: consider moving this outside of assert_energies2()
        system = deepcopy(system)
        for index in range(system.getNumForces()):
            force = system.getForce(index)
            force.setForceGroup(index) 
        integrator = openmm.VerletIntegrator(1.0*u.femtosecond)
        platform = openmm.Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)
        context.setPositions(positions)
        energy_components = list()
        for index in range(system.getNumForces()):
            force = system.getForce(index)
            forcename = force.__class__.__name__
            groups = 1<<index
            potential = beta * context.getState(getEnergy=True, groups=groups).getPotentialEnergy()
            energy_components.append((forcename, potential))
        del context, integrator
        return energy_components

    amber_energies = compute_potential_components(system_amber, parm_amber.positions)
    omm_energies = compute_potential_components(system_omm, parm_amber.positions)
    for amber, omm in zip(amber_energies, omm_energies):
        force_name = amber[0]
        print(force_name, amber[1], omm[1])

def assert_energies(prmtop, inpcrd, ffxml, system_name='unknown', tolerance=2.5e-5,
                    improper_tolerance=1e-2, units=u.kilojoules_per_mole, openmm_topology=None, openmm_positions=None):
    # AMBER
    parm_amber = parmed.load_file(prmtop, inpcrd)
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    amber_energies = parmed.openmm.energy_decomposition_system(parm_amber,
                     system_amber, nrg=units)

    # OpenMM-ffxml
    if isinstance(ffxml, str):
        ff = app.ForceField(ffxml)
    else:
        ff = app.ForceField(*ffxml)

    if openmm_positions is None:
        openmm_positions = parm_amber.positions

    if openmm_topology is not None:
        system_omm = ff.createSystem(openmm_topology)
        #parm_omm = parmed.openmm.load_topology(openmm_topology, system_omm,
        #                                       xyz=openmm_positions)
    else:
        system_omm = ff.createSystem(parm_amber.topology)
       # parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm,
       #                                        xyz=parm_amber.positions)
    #system_omm = parm_omm.createSystem(splitDihedrals=True)
    omm_energies = parmed.openmm.energy_decomposition_system(parm_amber,
                   system_omm, nrg=units, platform='Reference')
    print("amber energies", amber_energies)
    print("omm energies", omm_energies)
    # calc rel energies and assert
    energies = []
    rel_energies = []
    for i, j in zip(amber_energies, omm_energies):
        if i[0] != j[0]:
            raise Exception('Mismatch in energy tuples naming.')
        if abs(i[1]) > NEARLYZERO:
            rel_energies.append((i[0], abs((i[1]-j[1])/i[1])))
        else:
            if abs(j[1]) > NEARLYZERO:
                raise AssertionError('One of AMBER %s energies (%s) for %s is zero, '
                      'while the corresponding OpenMM energy is non-zero' %
                      (system_name, i[0], ffxml))
            rel_energies.append((i[0], 0))

    dihedrals_done = False
    for (i, amber_energy, openmm_energy) in zip(rel_energies, amber_energies, omm_energies):
        if i[0] != 'PeriodicTorsionForce':
            if i[1] > tolerance:
                raise AssertionError('%s relative energy error (%s, %f) outside of allowed tolerance (%f) for %s: AMBER %s OpenMM %s' %
                                     (system_name, i[0], i[1], tolerance, ffxml, amber_energy, openmm_energy))
        else:
            if not dihedrals_done:
                if i[1] > tolerance:
                    raise AssertionError('%s relative energy error (%s, %f) outside of allowed tolerance (%f) for %s: AMBER %s OpenMM %s' %
                                         (system_name, i[0], i[1], tolerance, ffxml, amber_energy, openmm_energy))
                dihedrals_done = True
            else: #impropers
                if i[1] > improper_tolerance:
                    raise AssertionError('%s relative energy error (%s-impropers, %f) outside of allowed tolerance (%f) for %s: AMBER %s OpenMM %s' %
                                         (system_name, i[0], i[1], improper_tolerance, ffxml, amber_energy, openmm_energy))

    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log['ffxml_name'] = ffxml
        amber_energies_log['test_system'] = system_name
        amber_energies_log['data_type'] = 'AMBER'
        amber_energies_log['units'] = units
        omm_energies_log['ffxml_name'] = ffxml
        omm_energies_log['test_system'] = system_name
        omm_energies_log['data_type'] = 'OpenMM'
        omm_energies_log['units'] = units
        rel_energies_log['ffxml_name'] = ffxml
        rel_energies_log['test_system'] = system_name
        rel_energies_log['data_type'] = 'abs((AMBER-OpenMM)/AMBER)'
        dihedrals_done = False
        for item in amber_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                amber_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                amber_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                amber_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in omm_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                omm_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in rel_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                rel_energies_log[item[0]] = item[1]

        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)

def validate_protein(ffxml_name, leaprc_name, united_atom=False):
    if verbose: print('Protein energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    ala3_top = tempfile.mkstemp()
    ala3_crd = tempfile.mkstemp()
    villin_top = tempfile.mkstemp()
    villin_crd = tempfile.mkstemp()
    leap_script_ala3_file = tempfile.mkstemp()
    leap_script_villin_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    if not united_atom:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])
    else:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3_ua.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin_ua.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])

    write_file(leap_script_ala3_file[0], leap_script_ala3_string)
    write_file(leap_script_villin_file[0], leap_script_villin_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_ala3_file[1], os.devnull))
    if os.path.getsize(ala3_top[1]) == 0 or os.path.getsize(ala3_crd[1]) == 0:
        raise LeapException(leap_script_ala3_file[1])
    os.system('tleap -f %s > %s' % (leap_script_villin_file[1], os.devnull))
    if os.path.getsize(villin_top[1]) == 0 or os.path.getsize(villin_crd[1]) == 0:
        raise LeapException(leap_script_villin_file[1])

    try:
        if verbose: print('Calculating and validating ala_ala_ala energies...')
        assert_energies(ala3_top[1], ala3_crd[1], ffxml_name,
                        system_name='protein-ala_ala_ala')
        if verbose: print('Ala_ala_ala energy validation successful!')

        if verbose: print('Calculating and validating villin headpiece energies...')
        assert_energies(villin_top[1], villin_crd[1], ffxml_name,
                        system_name='protein-villin headpiece')
        if verbose: print('Villin headpiece energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (ala3_top, ala3_crd, villin_top, villin_crd, leap_script_ala3_file,
                 leap_script_villin_file):
            os.unlink(f[1])
    if verbose: print('Protein energy validation for %s done!' % ffxml_name)

def validate_dna(ffxml_name, leaprc_name):
    if verbose: print('DNA energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    dna_top = tempfile.mkstemp()
    dna_crd = tempfile.mkstemp()
    leap_script_dna_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    leap_script_dna_string = """addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "DG" "DG5"  } { 1 "DG" "DG3"  }
{ 0 "DA" "DA5"  } { 1 "DA" "DA3"  }
{ 0 "DC" "DC5"  } { 1 "DC" "DC3"  }
{ 0 "DT" "DT5"  } { 1 "DT" "DT3"  }
}
x = loadPdb files/4rzn_dna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, dna_top[1], dna_crd[1])

    write_file(leap_script_dna_file[0], leap_script_dna_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_dna_file[1], os.devnull))
    if os.path.getsize(dna_top[1]) == 0 or os.path.getsize(dna_crd[1]) == 0:
        raise LeapException(leap_script_dna_file[1])

    try:
        if verbose: print('Calculating and validating DNA energies...')
        assert_energies(dna_top[1], dna_crd[1], ffxml_name,
                        system_name='nucleic-DNA')
        if verbose: print('DNA energy validation successful!')

    finally:
        if verbose: print('Deleting temp files...')
        for f in (dna_top, dna_crd, leap_script_dna_file):
            os.unlink(f[1])
    if verbose: print('DNA energy validation for %s done!' % ffxml_name)

def validate_rna(ffxml_name, leaprc_name):
    if verbose: print('RNA energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    rna_top = tempfile.mkstemp()
    rna_crd = tempfile.mkstemp()
    leap_script_rna_file = tempfile.mkstemp()
    leap_script_rna_file_alt = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')

    leap_script_rna_string = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "G5"  } { 1 "G" "G3"  } { "G" "G" }
{ 0 "A" "A5"  } { 1 "A" "A3"  } { "A" "A" }
{ 0 "C" "C5"  } { 1 "C" "C3"  } { "C" "C" }
{ 0 "U" "U5"  } { 1 "U" "U3"  } { "U" "U" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    leap_script_rna_string_alt = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "RG5"  } { 1 "G" "RG3"  } { "G" "RG" }
{ 0 "A" "RA5"  } { 1 "A" "RA3"  } { "A" "RA" }
{ 0 "C" "RC5"  } { 1 "C" "RC3"  } { "C" "RC" }
{ 0 "U" "RU5"  } { 1 "U" "RU3"  } { "U" "RU" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    write_file(leap_script_rna_file[0], leap_script_rna_string)
    write_file(leap_script_rna_file_alt[0], leap_script_rna_string_alt)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_rna_file[1], os.devnull))
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        # try alternative name mappings
        os.system('tleap -f %s > %s' % (leap_script_rna_file_alt[1], os.devnull))
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        raise LeapException(leap_script_rna_file_alt[1])

    try:
        if verbose: print('Calculating and validating RNA energies...')
        # improper testing turned off pending solution to problems
        assert_energies(rna_top[1], rna_crd[1], ffxml_name,
                        system_name='nucleic-RNA')
        if verbose: print('RNA energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (rna_top, rna_crd, leap_script_rna_file, leap_script_rna_file_alt):
            os.unlink(f[1])
    if verbose: print('RNA energy validation for %s done!' % ffxml_name)

def validate_gaff(ffxml_name, leaprc_name, gaff_dat_name):
    if verbose: print('GAFF energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    imatinib_top = tempfile.mkstemp()
    imatinib_crd = tempfile.mkstemp()
    leap_script_imatinib_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    leap_script_imatinib_string = """\
source %s
loadamberparams gaff/dat/%s
loadamberparams files/frcmod.imatinib
x = loadMol2 files/imatinib.mol2
saveAmberParm x %s %s
quit""" % (leaprc_name, gaff_dat_name, imatinib_top[1], imatinib_crd[1])
    write_file(leap_script_imatinib_file[0], leap_script_imatinib_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_imatinib_file[1], os.devnull))
    if os.path.getsize(imatinib_top[1]) == 0 or os.path.getsize(imatinib_crd[1]) == 0:
        raise LeapException(leap_script_imatinib_file[1])

    try:
        if verbose: print('Calculating and validating imatinib energies...')
        assert_energies(imatinib_top[1], imatinib_crd[1], (ffxml_name,
                        'files/imatinib.xml', 'files/imatinib_frcmod.xml'),
                         system_name='gaff-imatinib')
        if verbose: print('Imatinib energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (imatinib_top, imatinib_crd, leap_script_imatinib_file):
            os.unlink(f[1])
    if verbose: print('GAFF energy validation for %s done!' % ffxml_name)

def validate_phospho_protein(ffxml_name, leaprc_name,
                             supp_leaprc_name = 'oldff/leaprc.ff99SBildn',
                             supp_ffxml_name='ffxml/ff99SBildn.xml',
                             phospho="phospho10"):
    if '14' in leaprc_name:
        # Use AMBER14SB
        supp_leaprc_name = 'oldff/leaprc.ff14SB'
        supp_ffxml_name='ffxml/ff14SB.xml'
        phospho = 'phospho14'

    # this function assumes ffxml/ff14SB.xml already exists
    if verbose: print('Phosphorylated protein energy validation for %s' %
                      ffxml_name)
    for pdbname in glob.iglob(f'files/{phospho}/*.pdb'):
        if verbose: print('Now testing with pdb %s' % os.path.basename(pdbname))
        if verbose: print('Preparing temporary files for validation...')
        top = tempfile.mkstemp()
        crd = tempfile.mkstemp()
        leap_script_file = tempfile.mkstemp()

        if verbose: print('Preparing LeaP scripts...')
        leap_script_string = """source %s
source %s
x = loadPdb %s
saveAmberParm x %s %s
quit""" % (supp_leaprc_name, leaprc_name, pdbname, top[1], crd[1])

        write_file(leap_script_file[0], leap_script_string)

        if verbose: print('Running LEaP...')
        os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
        if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
            raise LeapException(leap_script_file[1])

        try:
            if verbose: print('Calculating and validating energies...')
            assert_energies(top[1], crd[1], (supp_ffxml_name, ffxml_name),
                            system_name='phospho_protein: %s'
                            % os.path.basename(pdbname))
            if verbose: print('Energy validation successful!')
        finally:
            if verbose: print('Deleting temp files...')
            for f in (top, crd, leap_script_file):
                os.unlink(f[1])
        if verbose: print('Phosphorylated protein energy validation for %s done!'
                          % ffxml_name)

def validate_water_ion(ffxml_name, source_recipe_files, solvent_name, recipe_name,
                       standard_ffxml=None):
    if verbose: print('Water and ions energy validation for %s' %
                      ffxml_name)
    if solvent_name == 'tip3p':
        HOH = 'TP3'
        solvent_frcmod = None
    elif solvent_name == 'tip4pew':
        HOH = 'T4E'
        solvent_frcmod = 'frcmod.tip4pew'
    elif solvent_name == 'spce':
        HOH = 'SPC'
        solvent_frcmod = 'frcmod.spce'
    elif solvent_name == 'tip3pfb':
        HOH = 'FB3'
        solvent_frcmod = 'frcmod.tip3pfb'
    elif solvent_name == 'tip4pfb':
        HOH = 'FB4'
        solvent_frcmod = 'frcmod.tip4pfb'
    pdb_name = 'files/water_ion/' + recipe_name + '.pdb'
    if verbose: print('Preparing temporary files for validation...')
    top = tempfile.mkstemp()
    crd = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()
    if verbose: print('Preparing LeaP scripts...')
    leap_script_string_part1 = """loadamberparams parm10.dat
loadamberparams %s
loadamberparams %s\n""" % (source_recipe_files[0], source_recipe_files[1])

    leap_script_string_part2 = """\nloadOff atomic_ions.lib
loadoff solvents.lib
HOH = %s
# for TIP4PEW
addPdbAtomMap {{ "M" "EPW" }}
x = loadPdb %s
saveAmberParm x %s %s
quit""" % (HOH, pdb_name, top[1], crd[1])

    if solvent_frcmod:
        leap_script_string = (leap_script_string_part1 + ('loadamberparams %s'
                               % solvent_frcmod) + leap_script_string_part2)
    else:
        leap_script_string = leap_script_string_part1 + leap_script_string_part2

    write_file(leap_script_file[0], leap_script_string)

    # this test does it's own energy assertion because of differences
    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
    if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
        raise LeapException(leap_script_file[1])
    try:
        if verbose: print('Calculating and validating energies...')
        pdb = app.PDBFile(pdb_name, extraParticleIdentifier='')
        if standard_ffxml is None:
            ff = app.ForceField(ffxml_name)
        else:
            ff = app.ForceField(ffxml_name, standard_ffxml)
        system_omm = ff.createSystem(pdb.topology)
        parm_omm = parmed.openmm.load_topology(pdb.topology, xyz=pdb.positions)
        parm_amber = parmed.load_file(top[1], crd[1])
        system_amber = parm_amber.createSystem()
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm,
                       system_omm, nrg=u.kilojoules_per_mole)
        for entry in omm_energies:
            if entry[0] == 'NonbondedForce':
                omm_nonbonded = entry[1]
        amber_energies = parmed.openmm.energy_decomposition_system(parm_amber,
                         system_amber, nrg=u.kilojoules_per_mole)
        for entry in amber_energies:
            if entry[0] == 'NonbondedForce':
                amber_nonbonded = entry[1]

        rel_nonbonded = abs((amber_nonbonded-omm_nonbonded) / amber_nonbonded)

        if rel_nonbonded > 1e-5:
            raise AssertionError('NonbondedForce Water and ions energy (%f) outside of '
                                 'allowed tolerance (%f) for %s:' % (rel_nonbonded, 1e-5, ffxml_name))
        if verbose: print('Energy validation successful!')

    finally:
        if verbose: print('Deleting temp files...')
        for f in (top, crd, leap_script_file):
            os.unlink(f[1])
    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log['ffxml_name'] = ffxml_name
        amber_energies_log['test_system'] = 'water_ion'
        amber_energies_log['data_type'] = 'AMBER'
        amber_energies_log['NonbondedForce'] = amber_nonbonded
        amber_energies_log['units'] = u.kilojoules_per_mole
        omm_energies_log['ffxml_name'] = ffxml_name
        omm_energies_log['test_system'] = 'water_ion'
        omm_energies_log['data_type'] = 'OpenMM'
        omm_energies_log['NonbondedForce'] = omm_nonbonded
        omm_energies_log['units'] = u.kilojoules_per_mole
        rel_energies_log['ffxml_name'] = ffxml_name
        rel_energies_log['test_system'] = 'water_ion'
        rel_energies_log['data_type'] = 'abs((AMBER-OpenMM)/AMBER)'
        rel_energies_log['NonbondedForce'] = rel_nonbonded
        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)
    if verbose: print('Water and ions energy validation for %s done!'
                      % ffxml_name)

def validate_impropers(ffxml_name, leaprc_name):
    if verbose: print('Impropers validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    top_villin = tempfile.mkstemp()
    crd_villin = tempfile.mkstemp()
    top_dna = tempfile.mkstemp()
    crd_dna = tempfile.mkstemp()
    top_rna = tempfile.mkstemp()
    crd_rna = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    leap_script_string = """source %s
x = loadPdb files/villin.pdb
y = loadPdb files/4rzn_dna.pdb
z = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
saveAmberParm y %s %s
saveAmberParm z %s %s
quit""" % (leaprc_name, top_villin[1], crd_villin[1], top_dna[1], crd_dna[1],
           top_rna[1], crd_rna[1])
    write_file(leap_script_file[0], leap_script_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
    if os.path.getsize(top_villin[1]) == 0 or os.path.getsize(crd_villin[1]) == 0:
        raise LeapException(leap_script_file[1])
    if os.path.getsize(top_dna[1]) == 0 or os.path.getsize(crd_dna[1]) == 0:
        raise LeapException(leap_script_file[1])
    if os.path.getsize(top_rna[1]) == 0 or os.path.getsize(crd_rna[1]) == 0:
        raise LeapException(leap_script_file[1])

    # load into parmed
    parm_amber_villin = parmed.load_file(top_villin[1])
    parm_amber_dna = parmed.load_file(top_dna[1])
    parm_amber_rna = parmed.load_file(top_rna[1])

    # OpenMM
    ff = app.ForceField(ffxml_name)
    sys_omm_villin = ff.createSystem(parm_amber_villin.topology)
    sys_omm_dna = ff.createSystem(parm_amber_dna.topology)
    sys_omm_rna = ff.createSystem(parm_amber_rna.topology)
    parm_omm_villin = parmed.openmm.load_topology(parm_amber_villin.topology,
                                                  sys_omm_villin)
    parm_omm_dna = parmed.openmm.load_topology(parm_amber_dna.topology,
                                               sys_omm_dna)
    parm_omm_rna = parmed.openmm.load_topology(parm_amber_rna.topology,
                                               sys_omm_rna)

    # prepare sets of idxs
    set_amber_villin = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_amber_villin.dihedrals if dih.improper])
    set_amber_dna = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_amber_dna.dihedrals if dih.improper])
    set_amber_rna = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_amber_rna.dihedrals if dih.improper])
    set_omm_villin = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_omm_villin.dihedrals if dih.improper])
    set_omm_dna = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_omm_dna.dihedrals if dih.improper])
    set_omm_rna = set([(dih.atom1.idx, dih.atom2.idx, dih.atom3.idx,
        dih.atom4.idx) for dih in parm_omm_rna.dihedrals if dih.improper])

    try:
        if (set_amber_villin - set_omm_villin != set() or
       set_omm_villin - set_amber_villin != set()):
            raise AssertionError("""Impropers validation fail for %s (villin)
                                    set_amber - set_omm: %s
                                    set_omm - set_amber: %s""" % (ffxml_name,
                                    set_amber_villin-set_omm_villin,
                                    set_omm_villin-set_amber_villin))
        if (set_amber_dna - set_omm_dna != set() or
       set_omm_dna - set_amber_dna != set()):
            raise AssertionError("""Impropers validation fail for %s (DNA)
                                    set_amber - set_omm: %s
                                    set_omm - set_amber: %s""" % (ffxml_name,
                                    set_amber_dna-set_omm_dna,
                                    set_omm_dna-set_amber_dna))
        if (set_amber_rna - set_omm_rna != set() or
       set_omm_rna - set_amber_rna != set()):
            raise AssertionError("""Impropers validation fail for %s (RNA)
                                    set_amber - set_omm: %s
                                    set_omm - set_amber: %s""" % (ffxml_name,
                                    set_amber_rna-set_omm_rna,
                                    set_omm_rna-set_amber_rna))
    finally:
        if verbose: print('Deleting temp files...')
        for f in (top_villin, crd_villin, top_dna, crd_dna, top_rna, crd_rna,
                  leap_script_file):
            os.unlink(f[1])
    if verbose: print('Improper validation for %s done!' % ffxml_name)

def validate_lipids(ffxml_name, leaprc_name):
    if verbose: print('Lipids energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    lipids_top = tempfile.mkstemp()
    lipids_crd = tempfile.mkstemp()
    leap_script_lipids_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    leap_script_lipids_string = """source %s
x = loadPdb files/POPC-nowater-amber.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, lipids_top[1], lipids_crd[1])
    write_file(leap_script_lipids_file[0], leap_script_lipids_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_lipids_file[1], os.devnull))
    if os.path.getsize(lipids_top[1]) == 0 or os.path.getsize(lipids_crd[1]) == 0:
        raise LeapException(leap_script_lipids_file[1])

    try:
        if verbose: print('Calculating and validating lipids energies...')
        assert_energies(lipids_top[1], lipids_crd[1], ffxml_name,
                        system_name='lipids')
        if verbose: print('Lipids energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (lipids_top, lipids_crd, leap_script_lipids_file):
            os.unlink(f[1])
    if verbose: print('Lipids energy validation for %s done!' % ffxml_name)

def validate_merged_lipids(ffxml_name, leaprc_name):
    if verbose: print('Lipids (merged) energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    lipids_top = tempfile.mkstemp()
    lipids_crd = tempfile.mkstemp()
    leap_script_lipids_file = tempfile.mkstemp()
    pdbfile = app.PDBFile('files/POPC-nowater-charmm.pdb')

    if verbose: print('Preparing LeaP scripts...')
    leap_script_lipids_string = """source %s
x = loadPdb files/POPC-nowater-amber.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, lipids_top[1], lipids_crd[1])
    write_file(leap_script_lipids_file[0], leap_script_lipids_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_lipids_file[1], os.devnull))
    if os.path.getsize(lipids_top[1]) == 0 or os.path.getsize(lipids_crd[1]) == 0:
        raise LeapException(leap_script_lipids_file[1])

    try:
        if verbose: print('Calculating and validating lipids energies...')
        assert_energies(lipids_top[1], lipids_crd[1], ffxml_name,
                        system_name='lipids',
                        openmm_topology=pdbfile.topology, openmm_positions=pdbfile.positions)
        if verbose: print('Lipids energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (lipids_top, lipids_crd, leap_script_lipids_file):
            os.unlink(f[1])
    if verbose: print('Lipids energy validation for %s done!' % ffxml_name)

def validate_glyco_protein(ffxml_name, leaprc_name,
                            supp_leaprc_name = 'oldff/leaprc.ff14SB',
                             supp_ffxml_name='ffxml/ff14SB.xml'):

    # this function assumes ffxml/ff14SB.xml already exists
    if verbose: print('Glycosylated protein energy validation for %s' %
                      ffxml_name)
    for pdbname in glob.iglob(f'files/glycan/*.pdb'):
        if verbose: print('Now testing with pdb %s' % os.path.basename(pdbname))
        if verbose: print('Preparing temporary files for validation...')
        top = tempfile.mkstemp()
        crd = tempfile.mkstemp()
        leap_script_file = tempfile.mkstemp()

        if verbose: print('Preparing LeaP scripts...')
        leap_script_string = """source %s
source %s
mol = loadPdb %s
#  Add disulphide bridges

bond mol.336.SG mol.361.SG
bond mol.379.SG mol.432.SG
bond mol.391.SG mol.525.SG
bond mol.480.SG mol.488.SG

##  N343 glycans

bond mol.343.ND2 mol.528.C1 #  Bond N343 to GlcNAc
bond mol.528.O6 mol.537.C1
bond mol.528.O4 mol.529.C1
bond mol.529.O4 mol.530.C1
bond mol.530.O6 mol.531.C1
bond mol.531.O2 mol.532.C1
bond mol.532.O4 mol.533.C1
bond mol.530.O3 mol.534.C1
bond mol.534.O2 mol.535.C1
bond mol.535.O4 mol.536.C1

saveAmberParm mol %s %s
quit""" % (leaprc_name, supp_leaprc_name,  pdbname, top[1], crd[1])

        write_file(leap_script_file[0], leap_script_string)

        if verbose: print('Running LEaP...')
        os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
        if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
            raise LeapException(leap_script_file[1])

        try:
            if verbose: print('Calculating and validating energies...')
            assert_energies2(top[1], crd[1], (supp_ffxml_name, ffxml_name),
                            system_name='glyco_protein: %s'
                            % os.path.basename(pdbname))
            if verbose: print('Energy validation successful!')
        finally:
            if verbose: print('Deleting temp files...')
            for f in (top, crd, leap_script_file):
                os.unlink(f[1])
        if verbose: print('Glycosylated protein energy validation for %s done!'
                          % ffxml_name)


class Logger():
    """
    Log energy discrepancies to a file.

    Parameters
    ----------
    log_filename : str
        Name of CSV file to write to

    """
    # logs testing energies into csv
    def __init__(self, log_filename=None):
        if log_filename:
            csvfile = open(log_filename, 'w')
            fieldnames = ['ffxml_name', 'data_type', 'test_system', 'units', 'HarmonicBondForce',
                      'HarmonicAngleForce', 'PeriodicTorsionForce_dihedrals',
                      'PeriodicTorsionForce_impropers', 'NonbondedForce']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            self.csvfile = csvfile
            self.writer = writer
        else:
            self.csvfile = None
            self.writer = None

    def close(self):
        if self.csvfile:
            self.csvfile.close()

    def log(self, energies):
        if self.writer:
            self.writer.writerow(energies)

if __name__ == '__main__':
    main()
