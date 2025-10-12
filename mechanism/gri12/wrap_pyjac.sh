#!/bin/bash

echo Preparing pyjac python wrapper for gri12.cti file

python3 -m pyjac --lang c --input gri12.cti 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
