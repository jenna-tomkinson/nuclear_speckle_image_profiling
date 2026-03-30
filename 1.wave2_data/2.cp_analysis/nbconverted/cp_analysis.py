#!/usr/bin/env python
# coding: utf-8

# # Perform nucleus segmentation and feature extraction
# 
# Using a CellProfiler pipeline, we will segment single nuclei and extract features across all channels.
# All data is outputted into a SQLite file.

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
run_name = "analysis"

# path to IC pipeline
path_to_pipeline = pathlib.Path("./pipeline/analysis.cppipe").resolve(strict=True)

# set main output dir for all plates
output_dir = pathlib.Path("./analysis_output")
output_dir.mkdir(exist_ok=True)

# directory where IC corrected images are located within folders
images_dir = pathlib.Path("../1.illumination_correction/IC_corrected_images").resolve(
    strict=True
)


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


# ## Run analysis pipeline on each plate in parallel
# 
# In this notebook, we do not run the cells to completion as we prefer to run the notebooks as nbconverted python files due to better stability.

# In[ ]:


cp_parallel.run_cellprofiler_parallel(
    plate_info_dictionary=plate_info_dictionary, run_name=run_name
)

