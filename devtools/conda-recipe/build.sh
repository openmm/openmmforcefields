#!/bin/bash

# Install the conversion tools package
#$PYTHON setup.py install

# Package the CHARMM forcefield
CHARMM_FFXML_DEST=""$PREFIX/share/openmm/forcefields/charmm"
mkdir -p $CHARMM_FFXML_DEST
cp -r charmm/ffxml/*.xml $CHARMM_FFXML_DEST

# Package the AMBER forcefield
