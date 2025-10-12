#!/bin/bash

echo Preparing pyjac python wrapper for KAUST_32s243r.inp file

python3 -m pyjac --lang c --input KAUST_32s243r.inp --thermo therm.dat -ls N2 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
