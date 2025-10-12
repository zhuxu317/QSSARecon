#!/bin/bash

echo Preparing pyjac python wrapper for Mei_39s256r.cti file

python3 -m pyjac --lang c --input Mei_39s256r.cti --last_species AR  2> err_lib.txt

python3 -m pyjac.pywrap --source_dir ./out/ --lang c 2> err_wrap.txt

  
