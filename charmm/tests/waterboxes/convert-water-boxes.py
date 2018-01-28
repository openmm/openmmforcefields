from parmed.charmm import CharmmParameterSet
from parmed.formats.psf import PSFFile
import parmed as pmd

type_map = {
    'O' : 'OT',
    'OH2' : 'OT',
    'H1' : 'HT',
    'H2' : 'HT',
    'LP' : 'LP',
    'OM' : 'LP',
    'LP1' : 'LP',
    'LP2' : 'LP',
}

for nsites in [3, 4, 5]:
    pdb_filename = 'waterbox-{}-site.pdb'.format(nsites)
    psf_filename = 'waterbox-{}-site.psf'.format(nsites)
    str_filename = 'toppar_tip{}p.str'.format(nsites)
    pdb = pmd.load_file(pdb_filename)
    for atom in pdb.atoms:
        atom.type = type_map[atom.name]
    # Delete bonds to lone pairs
    indices_to_delete = list()
    for (index, bond) in enumerate(list(pdb.bonds)):
        if (bond.atom1.type == 'LP') or (bond.atom2.type == 'LP'):
            indices_to_delete.append(index)
    for index in indices_to_delete[-1:0:-1]:
        pdb.bonds[index].delete()
        del pdb.bonds[index]

    params = CharmmParameterSet(str_filename)
    PSFFile.write(pdb, psf_filename)
