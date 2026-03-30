#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the main conda environment
conda activate nuclear_speckle_cp_env

# convert notebook to python file into the nbconverted folder
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

### Whole image qc was performed manually ###

# run python script to perform ic
python nbconverted/cp_ic.py
