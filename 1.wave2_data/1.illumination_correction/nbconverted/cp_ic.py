#!/usr/bin/env python
# coding: utf-8

# # Process TIFF files with illumination correction
# 
# Illumination correction (IC) is a step in the image analysis pipeline to remove any illumination errors that could be impacting the biology of the cells. 
# The function we generate in the CellProfiler pipeline are the "best" possible given the parameters available.
# It is run and corrected per image per channel.
# We will process the 16-bit TIFF files with IC and save as the same bit-depth.

# ## Import libraries

# In[1]:


import pathlib
import pprint

import sys

sys.path.append("../../utils")
import cp_parallel


# ## Set paths and variables

# In[2]:


# set the run type for the parallelization
run_name = "illum_correction"

# path to IC pipeline
path_to_pipeline = pathlib.Path("./pipeline/illumination_correction.cppipe").resolve(
    strict=True
)

# set main output dir for all plates
output_dir = pathlib.Path("./IC_corrected_images")
output_dir.mkdir(exist_ok=True)

# directory where images are located within folders
images_dir = pathlib.Path("../0.download_data/images").resolve(strict=True)


# ## Create dictionary with all info for each plate

# In[3]:


# create plate info dictionary with all parts of the CellProfiler CLI command to run in parallel
plate_info_dictionary = {
    "wave2_data": {
        "path_to_images": images_dir,
        "path_to_output": output_dir,
        "path_to_pipeline": path_to_pipeline,
    }
}

# view the dictionary to assess that all info is added correctly
pprint.pprint(plate_info_dictionary, indent=4)


# ## Run illumination correction pipeline on each plate in parallel
# 
# In this notebook, we do not run the cells to completion as we prefer to run the notebooks as nbconverted python files due to better stability.

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)

