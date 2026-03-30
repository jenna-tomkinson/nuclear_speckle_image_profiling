#!/usr/bin/env python
# coding: utf-8

# # Convert SQLite file output from CellProfiler into parquet file using Cytotable

# ## Import libraries

# In[1]:


import logging
import pathlib

import pandas as pd

# cytotable will merge objects from SQLite file into single cells and save as parquet file
from cytotable import convert, presets

# Set the logging level to a higher level to avoid outputting unnecessary errors from config file in convert function
logging.getLogger().setLevel(logging.ERROR)


# ## Set paths and variables

# In[2]:


# type of file output for CytoTable
dest_datatype = "parquet"

# set main output dir for all parquet files
output_dir = pathlib.Path("./data")
output_dir.mkdir(exist_ok=True)

# directory where SQLite files are located
sqlite_dir = pathlib.Path("../2.cp_analysis/analysis_output").resolve(strict=True)

# Set converted parquet dir
parquet_dir = pathlib.Path(f"{output_dir}/converted_profiles")
parquet_dir.mkdir(exist_ok=True)


# ## Run cytotable convert to output nuclei and image features separately for all plates

# In[3]:


# Construct output path for converted parquet file
output_path = pathlib.Path(f"{parquet_dir}/u2os_data_converted.parquet")

print("Starting conversion with CytoTable for dataset")

# merge single cells and output as parquet file
convert(
    source_path=str(sqlite_dir),
    dest_path=str(output_path),
    dest_datatype=dest_datatype,
    metadata=["image"],
    compartments=["nuclei"],
    identifying_columns=["ImageNumber"],
    page_keys={
        "image": "ImageNumber",
        "nuclei": "Nuclei_Number_Object_Number",
        "join": "Nuclei_Number_Object_Number",
    },
    joins="""
    SELECT
        *
    FROM
        read_parquet('per_image.parquet') as per_image
    INNER JOIN read_parquet('per_nuclei.parquet') AS per_nuclei ON
        per_nuclei.Metadata_ImageNumber = per_image.Metadata_ImageNumber
""",
    chunk_size=10000,
)

print("Conversion finished for dataset")


# ## Remove unwanted image + metadata columns and split the bulk and single-cell data from the main parquet file

# In[4]:


# path to unwanted image cols text file
unwanted_list_path = pathlib.Path("./unwanted_image_cols.txt")

# Load the list of columns to remove from the text file
with open(unwanted_list_path, "r") as file:
    columns_to_remove = [line.strip() for line in file]

# Read in the single parquet file
plate_df = pd.read_parquet(parquet_dir / "u2os_data_converted.parquet")

print("Starting to edit image and nuclei data frames")

# Drop the specified columns (ignore error if a column isn't there)
plate_df = plate_df.drop(columns=columns_to_remove, errors="ignore")

# Identify metadata columns for nuclei data frame
metadata_columns = [
    "Metadata_ImageNumber",
    "Image_Metadata_Plate",
    "Image_Metadata_Position",
    "Image_Metadata_Best_Z",
    "Image_Metadata_Well",
    "Image_Count_Nuclei",
    "Image_FileName_DAPI",
    "Image_PathName_DAPI",
]

# Create nuclei (single-cell) data frame
nuclei_df = plate_df[
    metadata_columns + [col for col in plate_df.columns if col.startswith("Nuclei_")]
]

# Create image (bulk) data frame
image_df = plate_df[
    ["Metadata_ImageNumber"]
    + [col for col in plate_df.columns if col.startswith("Image_")]
]

# Drop duplicate images since each image repeats across cells
image_df = image_df.drop_duplicates(subset="Metadata_ImageNumber")

# Save outputs (same directory as input parquet)
nuclei_df.to_parquet(parquet_dir / "u2os_per_nuclei.parquet", index=False)
image_df.to_parquet(parquet_dir / "u2os_per_image.parquet", index=False)

# Quick sanity prints
print("Shape of nuclei data frame", nuclei_df.shape)
print("Shape of image data frame", image_df.shape)

