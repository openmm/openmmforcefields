#!/usr/bin/env bash

set -e

python convert_amber_ions.py "${AMBERHOME}"
python convert_amber.py --input solvents.yaml --verbose
python convert_amber.py --input biopolymer.yaml --verbose
python convert_amber.py --input gaff.yaml --verbose
python convert_amber.py --input glycam/glycan.yaml --verbose
python convert_amber.py --combination-tests --verbose
