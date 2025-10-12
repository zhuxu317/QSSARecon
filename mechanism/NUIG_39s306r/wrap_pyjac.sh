#!/bin/bash

echo Preparing pyjac python wrapper for NUIG_39s306r.cti file

python3 -m pyjac --lang c --input NUIG_39s306r.cti -fopt -ls N2 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
