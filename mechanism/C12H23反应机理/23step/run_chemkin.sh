#!/bin/bash

# Set up Chemkin environment
export CHEMKIN_ROOT=/opt/ansys_inc/v192/reaction/chemkinpro.linuxx8664
export PATH=$CHEMKIN_ROOT/bin:$PATH
export CHEMKIN_DATA=$CHEMKIN_ROOT/data

# Clean up old files
rm -f chem.asc tran.asc cklink ckinterp.out ckinterp.log

# Run Chemkin preprocessing
echo "Running Chemkin preprocessing..."
echo "jeta-23steps.inp
therm-23steps.dat
tran-23steps.dat" | CKPreProcess

# Check if preprocessing was successful
if [ ! -f chem.asc ] || [ ! -f tran.asc ]; then
    echo "ERROR: Chemkin preprocessing failed"
    exit 1
fi

# Run PREMIX for flame speed calculation
echo "Running PREMIX flame speed calculation..."
CKReactorFreelyPropagatingFlame < flame_speed.inp > flame_speed.out

echo "Calculation complete. Check flame_speed.out for results."