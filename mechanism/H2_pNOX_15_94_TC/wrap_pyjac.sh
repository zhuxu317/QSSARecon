#!/bin/bash

echo Preparing pyjac python wrapper for H2_pNOX_15_94_TC.cti file

python3 -m pyjac --lang c --input H2_pNOX_15_94_TC.cti --last_species AR 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
