#!/bin/bash

# Set up Chemkin environment
export CHEMKIN_ROOT=/opt/ansys_inc/v192/reaction/chemkinpro.linuxx8664
export PATH=$CHEMKIN_ROOT/bin:$PATH
export CHEMKIN_DATA=$CHEMKIN_ROOT/data

# Create directories for both mechanisms
mkdir -p results_clean results_original

# Function to run flame speed calculation
run_flame_speed() {
    local mechanism=$1
    local workdir=$2
    local label=$3
    
    echo "=== Running flame speed calculation for $label ==="
    
    cd $workdir
    
    # Clean up
    rm -f chem.asc tran.asc cklink ckinterp.out ckinterp.log flame_speed.out
    
    # Copy necessary files
    cp ../$mechanism .
    cp ../therm-23steps.dat .
    cp ../tran-23steps.dat .
    
    # Run preprocessing
    echo "Running Chemkin preprocessing..."
    echo "$mechanism
therm-23steps.dat
tran-23steps.dat" | CKPreProcess
    
    # Check if preprocessing succeeded
    if [ ! -f chem.asc ] || [ ! -f tran.asc ]; then
        echo "ERROR: Chemkin preprocessing failed for $label"
        cd ..
        return 1
    fi
    
    # Create flame speed input
    cat > flame_speed.inp <<EOF
FLAME SPEED CALCULATION FOR $label MECHANISM

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
    
    # Run flame speed calculation
    echo "Running flame speed calculation..."
    CKReactorFreelyPropagatingFlame < flame_speed.inp > flame_speed.out 2> error.log
    
    # Extract results
    if [ -f flame_speed.out ]; then
        echo "=== Results for $label ==="
        grep -i "flame.*speed\|velocity\|cm/s\|m/s" flame_speed.out | head -10
        tail -20 flame_speed.out
    else
        echo "ERROR: No output file generated for $label"
    fi
    
    cd ..
}

# Run calculations for both mechanisms
echo "Starting laminar flame speed calculations..."
echo "Working directory: $(pwd)"
echo ""

run_flame_speed "jeta_clean.inp" "results_clean" "CLEAN_MECHANISM"
echo ""
run_flame_speed "jeta-23steps.inp" "results_original" "ORIGINAL_MECHANISM"

echo ""
echo "=== Calculation Summary ==="
echo "Results saved in:"
echo "- results_clean/  (jeta_clean.inp)"
echo "- results_original/  (jeta-23steps.inp)"

# Create summary file
cat > results_summary.txt <<EOF
Laminar Flame Speed Calculation Summary
=====================================

Mechanisms tested:
1. jeta_clean.inp - Simplified Chemkin-compatible format
2. jeta-23steps.inp - Original mechanism format

Calculation conditions:
- Pressure: 1.0 atm
- Temperature: 300.0 K
- Fuel: C12H23 (φ = 1.0)
- Oxidizer: Air (O2:N2 = 3.5:13.16)

Results:
EOF

echo "Analysis complete. Check results_summary.txt for summary."