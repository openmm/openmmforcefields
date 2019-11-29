import xml.etree.ElementTree as etree
from copy import deepcopy

# Define file paths

charmm_ffxml_filename = '../charmm/ffxml/charmm36_nowaters.xml' # input CHARMM ffxml file containing lipid definitions
amber_unmerged_lipid_ffxml_filename = 'ffxml/lipid17.xml' # input AMBER ffxml file for unmerged lipids
charmmlipid2amber_filename = 'files/charmmlipid2amber.csv' # AMBER to CHARMM lipid residue conversion table
amber_merged_lipid_ffxml_filename = 'ffxml/lipid17_merged.xml' # output AMBER ffxml file for merged lipids

# Read the input files.

charmmff = etree.parse(charmm_ffxml_filename)
amberff = etree.parse(amber_unmerged_lipid_ffxml_filename)
charmmResidues = charmmff.getroot().find('Residues').findall('Residue')
amberResidues = amberff.getroot().find('Residues').findall('Residue')
amberResMap = {}
for res in amberResidues:
    atoms = dict((atom.attrib['name'], atom) for atom in res.findall('Atom'))
    amberResMap[res.attrib['name']] = atoms
translations = {}
with open(charmmlipid2amber_filename) as input:
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
    return newres

for residue in charmmResidues:
    copy = translateResidue(residue)
    if copy is not None:
        parentNode.append(copy)

# Write merged lipid ffxml file

amberff.write(amber_merged_lipid_ffxml_filename)
