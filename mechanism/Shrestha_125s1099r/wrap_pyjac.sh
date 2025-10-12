#!/bin/bash

echo Preparing pyjac python wrapper for Shrestha_125s1099r.cti file

python3 -m pyjac --lang c --input Shrestha_125s1099r.cti 2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
