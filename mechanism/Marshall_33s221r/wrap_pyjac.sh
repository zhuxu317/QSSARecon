#!/bin/bash


# zhuxu21@rserver /d/Z/C/N/m/Marshall_33s221r> python3 -m pyjac --help
# usage: __main__.py [-h] -l {c,cuda,fortran,matlab} -i INPUT [-t THERMO]
#                    [-ic INITIAL_CONDITIONS] [-co] [-nosmem] [-pshare]
#                    [-nb NUM_BLOCKS] [-nt NUM_THREADS] [-mt MULTI_THREAD]
#                    [-fopt] [-b BUILD_PATH] [-ls LAST_SPECIES] [-ad] [-sj]

# pyJac: Generates source code for analytical chemical Jacobians.

# optional arguments:
#   -h, --help            show this help message and exit
#   -l {c,cuda,fortran,matlab}, --lang {c,cuda,fortran,matlab}
#                         Programming language for output source files.
#   -i INPUT, --input INPUT
#                         Input mechanism filename (e.g., mech.dat).
#   -t THERMO, --thermo THERMO
#                         Thermodynamic database filename (e.g., therm.dat), or
#                         nothing if in mechanism.
#   -ic INITIAL_CONDITIONS, --initial-conditions INITIAL_CONDITIONS
#                         A comma separated list of initial initial conditions
#                         to set in the set_same_initial_conditions method.
#                         Expected Form: T,P,Species1=...,Species2=...,...
#                         Temperature in K Pressure in Atm Species in moles
#   -co, --cache-optimizer
#                         Attempt to optimize cache store/loading via use of a
#                         greedy selection algorithm. (Experimental)
#   -nosmem, --no-shared-memory
#                         Use this option to turn off attempted shared memory
#                         acceleration for CUDA.
#   -pshare, --prefer-shared
#                         Use this option to allocate more space for shared
#                         memory than the L1 cache for CUDA (not recommended).
#   -nb NUM_BLOCKS, --num-blocks NUM_BLOCKS
#                         The target number of blocks / sm for CUDA.
#   -nt NUM_THREADS, --num-threads NUM_THREADS
#                         The target number of threads / block for CUDA.
#   -mt MULTI_THREAD, --multi-threaded MULTI_THREAD
#                         The number of threads to use during the optimization
#                         process.
#   -fopt, --force-optimize
#                         Use this option to force a reoptimization of the
#                         mechanism (usually only happens when generating for a
#                         different mechanism).
#   -b BUILD_PATH, --build_path BUILD_PATH
#                         The folder to generate the Jacobian and rate
#                         subroutines in.
#   -ls LAST_SPECIES, --last_species LAST_SPECIES
#                         The name of the species to set as the last in the
#                         mechanism. If not specifed, defaults to the first of
#                         N2, AR, and HE in the mechanism.
#   -ad, --auto_diff      Use this option to generate file for use with the
#                         Adept autodifferentiation library.
#   -sj, --skip_jac       If specified, this option turns off Jacobian
#                         generation (only rate subs are generated)



echo Preparing pyjac python wrapper for Marshall_33s221r.cti file

python3 -m pyjac --lang c --input Marshall_33s221r.cti --last_species AR  2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
