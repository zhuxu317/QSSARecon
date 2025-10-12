#!/bin/bash

echo Preparing pyjac python wrapper for Otomo_32s213r.cti file

python3 -m pyjac --lang c --input Otomo_32s213r.cti -fopt -LAST_SPECIES N2 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
