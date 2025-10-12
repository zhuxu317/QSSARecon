#!/bin/bash

echo Preparing pyjac python wrapper for Li_demo.cti file

python3 -m pyjac --lang c --input Li_demo.cti  2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt
