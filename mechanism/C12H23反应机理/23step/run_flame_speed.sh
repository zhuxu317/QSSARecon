#!/bin/bash

# Set up Chemkin environment
export CHEMKIN_ROOT=/opt/ansys_inc/v192/reaction/chemkinpro.linuxx8664
export PATH=$CHEMKIN_ROOT/bin:$PATH
export CHEMKIN_DATA=$CHEMKIN_ROOT/data

# Create working directory
mkdir -p chemkin_work
cd chemkin_work

# Copy mechanism files
cp ../jeta-23steps.inp .
cp ../therm-23steps.dat .
cp ../tran-23steps.dat .

# Clean up old files
rm -f chem.asc tran.asc cklink ckinterp.out ckinterp.log

echo "Starting Chemkin flame speed calculation for jeta-23steps mechanism..."
echo ""
echo "Step 1: Running Chemkin preprocessing..."

# Create preprocessing input file
cat > ck_preprocess_input.txt <<EOF
jeta-23steps.inp
therm-23steps.dat
tran-23steps.dat
EOF

# Run preprocessing
CKPreProcess < ck_preprocess_input.txt

# Check if preprocessing succeeded
if [ ! -f chem.asc ] || [ ! -f tran.asc ]; then
    echo "ERROR: Chemkin preprocessing failed"
    echo "Checking mechanism file format..."
    head -20 ../jeta-23steps.inp
    exit 1
fi

echo ""
echo "Step 2: Creating flame speed calculation input..."

# Create flame speed input file
cat > flame_speed.inp <<EOF
FLAME SPEED CALCULATION FOR JETA-23STEPS MECHANISM

PROBLEM TYPE
  PREMIXED FLAMES
  FREE FLAME
  FLAME SPEED

PRESSURE
  1.0 ATM

TEMPERATURE
  300.0 K

FUEL
  C12H23          1.0

OXIDIZER
  O2              3.5
  N2              13.16

EQUIVALENCE RATIO
  1.0

TRANSPORT
  MIXTURE-AVERAGED

GRID PARAMETERS
  ADAPT 1
  GRAD 0.1
  CURV 0.2

SOLUTION PARAMETERS
  MAX GRID POINTS 1000
  MAX SUBITERATIONS 20
  ABSOLUTE TOLERANCE 1.0E-8
  RELATIVE TOLERANCE 1.0E-4

END
EOF

echo ""
echo "Step 3: Running flame speed calculation..."

# Run flame speed calculation
CKReactorFreelyPropagatingFlame < flame_speed.inp > flame_speed.out

echo ""
echo "Calculation complete!"
echo "Results saved in: chemkin_work/flame_speed.out"
echo ""
echo "Looking for flame speed result..."
grep -i "flame.*speed\|velocity" flame_speed.out || echo "Search flame_speed.out manually for results"