#!/bin/bash -x

# This script is for running the preparation workflow to create a pdb file that is compatible for parsing by openmm
# It reads in a pdb provided in the supplementary information from 10.1021/ja512751q always of the form:
# {target}_protein.pdb and runs pdb4amber and pdbfixer to generate a new pdb file with form {target}_protein_fixed.pdb


if [ $1 == "-h" ]; then
    echo "
        This script assumes that the first argument is an unprocessed pdb file of the form
        {target}_protein.pdb. It also assumes pdb4amber and pdbfixer is accessible from your \$PATH"
    exit
fi

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 {input {target}_protein.pdb}"
    echo "Error: This script assumes that the first argument is an unprocessed pdb file of the form
        {target}_protein.pdb. It also assumes pdb4amber and pdbfixer is accessible from your \$PATH"
    exit
fi

if [ -z $AMBERHOME ]; then
  echo "Error: AMBERHOME is not set."
  echo "Please set AMBERHOME in your environment"
fi

receptor=$1
basename=`echo $receptor | cut -d'.' -f1`

grep -v HET $receptor > ${basename}_nohet.pdb
pdb4amber -i ${basename}_nohet.pdb -o ${basename}_pdb4amber.pdb
pdbfixer ${basename}_pdb4amber.pdb --output ${basename}_fixed.pdb --keep-heterogens water --verbose

exit