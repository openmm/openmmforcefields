#!/bin/bash
# Script to compute energies inside docker container

# Change to mount directory
cd /mnt

# Add charmm to path
export PATH="/charmm/c40b1_gnu/exec/gnu:$PATH"

# Run comparison
charmm < dhfr.inp
