#!/usr/bin/env bash

set -e

for input_file in files/*.yaml; do
    python convert_charmm.py --input "${input_file}" --verbose
done
