from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed import openmm
import glob
import yaml
from collections import OrderedDict
import hashlib
import os
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
import argparse
import csv
import logging
import warnings

def main():
    global verbose
    global no_log
    global logger
    # Set up parser
    parser = argparse.ArgumentParser(description='CHARMM --> OpenMM forcefield conversion script')
    parser.add_argument('--input', '-i', default='files/waters.yaml',
                        help='path of the input file. Default: "files/waters.yaml"')
    parser.add_argument('--output-dir', '-od', help='path of the output directory. '
                        'Default: "ffxml/" for yaml, "./" for leaprc', default='ffxml/')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='turns verbosity on')
    parser.add_argument('--no-log', action='store_true',
                        help='turns logging of energies to log.csv off')
    args = parser.parse_args()
    verbose = args.verbose
    no_log = args.no_log

    if not no_log: logger = Logger('log.csv')
    convert_yaml(args.input, ffxml_dir=args.output_dir)
    if not no_log: logger.close()

def convert_yaml(yaml_filename, ffxml_dir):
    # Read YAML
    data = yaml.safe_load(open(yaml_filename, 'r'))
    source_pack = data[0]['sourcePackage']
    source_pack_ver = data[0]['sourcePackageVersion']

    for entry in data[1:]:
        ffxml_filename = os.path.join(ffxml_dir, entry['Destination'])
        print('Generating %s' % ffxml_filename)
        charmm_references = entry['References']
        source_files = entry['Source']

        # files that should be excluded from conversion.
        exclude_files = set()
        if ('exclude' in source_files) and (source_files['exclude'] is not None):
            exclude_files = set(source_files['exclude'])

        # charmm36 main top and par files
        charmm_files = list()
        if ('include' in source_files) and (source_files['include'] is not None):
            charmm_files = source_files['include']

        # add stream files
        for files in source_files['stream']:
            charmm_files.extend(glob.glob(files))

        # exclude files from conversion
        charmm_files = set(charmm_files) - exclude_files

        provenance = OrderedDict()
        source = provenance['Source'] = []
        for fi in charmm_files:
            source.append(OrderedDict())
            source[-1]['Source'] = fi
            md5 = hashlib.md5()
            with open(fi, 'rb') as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]['md5hash'] = md5
            source[-1]['sourcePackage'] = source_pack
            source[-1]['sourcePackageVersion'] = source_pack_ver

        references = provenance['Reference'] = []
        for ff in charmm_references:
            for cite in charmm_references[ff]:
                references.append(OrderedDict())
                if type(cite) is dict:
                    for key in cite.keys():
                        citation = cite[key]
                        references[-1]['Reference'] = citation
                        references[-1]['forcefield'] = ff
                        references[-1]['type'] = key
                else:
                    citation = cite
                    references[-1]['Reference'] = citation
                    references[-1]['forcefield'] = ff

        #generate recommended combination for charmm36
        if verbose: print('Loading CHARMM parameter sets %s...' % charmm_files)
        params = CharmmParameterSet(*charmm_files)

        if verbose: print('Converting parameters to OpenMM...')
        params_omm = openmm.OpenMMParameterSet.from_parameterset(params)

        # Set override level
        if 'Override' in entry:
            override_level = int(entry['Override'])
            if verbose: print('Setting residues and patches to override level %d...' % override_level)
            for name, residue in params_omm.residues.items():
                residue.override_level = override_level
            for name, patch in params_omm.patches.items():
                patch.override_level = override_level

        if verbose: print('Writing parameter set and compatible patches. This may take several minutes...')
        params_omm.write(ffxml_filename, provenance=provenance)

        # Try reading the forcefield back in to make sure it is valid
        if verbose: print('Verifying ffxml file integrity...')
        if 'TestInclude' in entry:
            ffxml_include = entry['TestInclude']
            forcefield = app.ForceField(ffxml_filename, *ffxml_include)
        else:
            forcefield = app.ForceField(ffxml_filename)
        if verbose: print('Verified.')

class Logger():
    # logs testing energies into csv
    def __init__(self, log_file):
        csvfile = open(log_file, 'w')
        fieldnames = ['ffxml_name', 'data_type', 'test_system', 'units', 'HarmonicBondForce',
                      'HarmonicAngleForce', 'PeriodicTorsionForce_dihedrals',
                      'PeriodicTorsionForce_impropers', 'NonbondedForce']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        self.csvfile = csvfile
        self.writer = writer

    def close(self):
        self.csvfile.close()

    def log(self, energies):
        self.writer.writerow(energies)

if __name__ == '__main__':
    main()
