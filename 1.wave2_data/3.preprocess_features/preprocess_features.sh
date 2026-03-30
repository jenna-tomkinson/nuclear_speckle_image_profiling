#!/bin/bash

# initialize the correct shell for your machine to allow conda to work (see README for note on shell names)
conda init bash
# activate the main conda environment
conda activate python_nuclear_speckle_env

# convert notebook to python file into the nbconverted folder
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# run python script to preprocess features
python nbconverted/0.convert_cytotable.py
python nbconverted/1.single_cell_qc.py
python nbconverted/2.single_cell_processing.py
