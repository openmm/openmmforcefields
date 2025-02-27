#!/usr/bin/env bash

set -e

python test_charmm.py tests/*/*.yaml --charmm --openmm-charmm --openmm-ffxml --openmm-ffxml-fix-impropers --perturb-replicates 5 --perturb-seed 291700478 "$@"
