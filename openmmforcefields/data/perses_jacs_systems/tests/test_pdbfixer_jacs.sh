#!/bin/bash -x

# This script is for testing the preparation workflow to create a pdb file that is compatible for parsing by openmm
# It reads in a pdb provided in the supplementary information from 10.1021/ja512751q always of the form:
# {target}_protein.pdb and runs pdb4amber and pdbfixer to generate a new pdb file with form {target}_protein_fixed.pdb
# for all targets in the source directory


if [ $1 == "-h" ]; then
    echo "
        This script is for testing unprocessed pdb files of the form
        {target}_protein.pdb. It also assumes pdb4amber and pdbfixer is accessible from your \$PATH"
    exit
fi

if [ -z $AMBERHOME ]; then
  echo "Error: AMBERHOME is not set."
  echo "Please set AMBERHOME in your environment"
fi

top=`pwd`
for protein_path in ../*/*_protein.pdb
do
    cd `dirname $protein_path`
    protein=`basename $protein_path`
    ../run_pdbfixer.sh $protein
    python $top/test_outputs.py -p *_fixed.pdb
    cd $top
done

echo "Testing finished, examine output to determine if any errors are present"
exit