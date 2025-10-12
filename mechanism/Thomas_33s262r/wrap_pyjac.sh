#!/bin/bash

echo Preparing pyjac python wrapper for Thomas_33s262r.cti file

python3 -m pyjac --lang c --input Thomas_33s262r.cti -fopt -ls N2 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt