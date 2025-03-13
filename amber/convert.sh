#!/usr/bin/env bash

set -e

python convert_amber_ions.py "${AMBERHOME}"
python convert_amber.py --input solvents.yaml --verbose
python convert_amber.py --input biopolymer.yaml --verbose
python convert_amber.py --input gaff.yaml --verbose

# GLYCAM automated conversion is broken (incorrect external bonds get
# generated): see https://github.com/openmm/openmmforcefields/pull/180.
# If a new version of the force field is released with AmberTools that we want
# to support, this will have to be corrected manually again, or a way to check
# it automatically can be devised.
# python convert_amber.py --input glycam/glycan.yaml --verbose
