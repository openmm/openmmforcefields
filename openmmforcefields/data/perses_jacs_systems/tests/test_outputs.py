import sys
from argparse import ArgumentParser
from simtk.openmm import app


if __name__ == "__main__":
    parser = ArgumentParser(epilog='''This program reads a pdb file and tests whether the file that is compatible 
                                      with the pdb parsing capabilities in openmm. The first argument should be the 
                                      pdb file to be used for testing.''')
    group = parser.add_argument_group("Required Arguments")
    group.add_argument('-p', '--protein', default=None, dest='protein', help='Protein pdb file to test compatibility'
                                                                             'with openmm', required=True)
    args = parser.parse_args()
    pdb = app.PDBFile(args.protein)
    modeller = app.Modeller(pdb.topology, pdb.positions)
    hs = [atom for atom in modeller.topology.atoms() if atom.element.symbol in ["H"]]
    modeller.delete(hs)
    forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
    try:
        modeller.addHydrogens(forcefield=forcefield)
    except ValueError as error:
        print(error)
        sys.exit("File %s was not properly prepared and cannot be parsed by openmm" % args.protein)

    sys.exit("Test passed! Output %s was prepared properly and can be parsed by openmm" % args.protein)