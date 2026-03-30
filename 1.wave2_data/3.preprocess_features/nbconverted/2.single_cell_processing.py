#!/usr/bin/env python
# coding: utf-8

# # Perform data preprocessing with pycytominer on single cell features
# 
# Note: Single cell is represented by only nuclei compartment features which was used to extract features across all three channels.

# ## Import libraries

# In[1]:


import pathlib

import pandas as pd
from pycytominer import normalize, feature_select


# ## Set paths and variables

# In[2]:


# Path to dir with converted profiles per plate (each plate as a folder)
converted_dir = pathlib.Path("./data/converted_profiles")

# Path to dir with `qc.parquet` with metadata for which nuclei failed qc
qc_dir = pathlib.Path("./data/single_cell_qc")

# path for plate map directory
platemap_dir = pathlib.Path("../0.download_data/metadata/platemaps")

# Output dir for the files to be saved to
output_dir = pathlib.Path("./data/single_cell_profiles")
output_dir.mkdir(exist_ok=True, parents=True)

# operations to perform for feature selection
feature_select_ops = [
    "variance_threshold",
    "correlation_threshold",
    "blocklist",
]


# ## Perform preprocessing on single cell features

# In[3]:


profile_path = pathlib.Path(f"{converted_dir}/u2os_per_nuclei.parquet")

print("Performing pycytominer pipeline for dataset")

output_normalized_file = pathlib.Path(
    f"{output_dir}/u2os_per_nuclei_sc_normalized.parquet"
)
output_feature_select_file = pathlib.Path(
    f"{output_dir}/u2os_per_nuclei_sc_feature_selected.parquet"
)

# Load data
profile_df = pd.read_parquet(profile_path)

qc_file = qc_dir / "qc.parquet"
qc_df = pd.read_parquet(qc_file)

# ---- Define join columns ----
join_cols = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Position",
    "Image_Metadata_Best_Z",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
]

# ---- Select all Metadata_cqc columns ----
qc_cqc_cols = [col for col in qc_df.columns if col.startswith("Metadata_cqc")]

qc_subset = qc_df[join_cols + qc_cqc_cols]

# ---- Merge QC into profile ----
profile_df = profile_df.merge(qc_subset, on=join_cols, how="left")

# ---- Drop rows where any QC column indicates failure (True) ----
qc_fail_cols = [col for col in qc_subset.columns if col.startswith("Metadata_cqc")]
annotated_df = profile_df[~profile_df[qc_fail_cols].any(axis=1)]
print(
    f"After QC filtering, {len(annotated_df)} nuclei remain out of {len(profile_df)} total nuclei."
)

# Rename columns using the rename() function
column_name_mapping = {
    "Image_Metadata_Plate": "Metadata_Plate",
    "Image_Metadata_Well": "Metadata_Well",
    "Image_Metadata_Position": "Metadata_Position",
    "Image_Metadata_Best_Z": "Metadata_Best_Z",
    "Image_Count_Nuclei": "Metadata_Nuclei_Site_Count",
    "Nuclei_Location_Center_X": "Metadata_Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y": "Metadata_Nuclei_Location_Center_Y",
    "Nuclei_AreaShape_BoundingBoxMaximum_X": "Metadata_Nuclei_AreaShape_BoundingBoxMaximum_X",
    "Nuclei_AreaShape_BoundingBoxMaximum_Y": "Metadata_Nuclei_AreaShape_BoundingBoxMaximum_Y",
    "Nuclei_AreaShape_BoundingBoxMinimum_X": "Metadata_Nuclei_AreaShape_BoundingBoxMinimum_X",
    "Nuclei_AreaShape_BoundingBoxMinimum_Y": "Metadata_Nuclei_AreaShape_BoundingBoxMinimum_Y",
}

annotated_df.rename(columns=column_name_mapping, inplace=True)

# Step 2: Normalization
normalized_df = normalize(
    profiles=annotated_df,
    method="mad_robustize",
    output_file=output_normalized_file,
    output_type="parquet",
)

# Step 3: Feature selection
feature_select(
    output_normalized_file,
    operation=feature_select_ops,
    output_file=output_feature_select_file,
    na_cutoff=0,
    output_type="parquet",
    blocklist_file="./blocklist_features.txt",
)
print(
    "Annotation, normalization, and feature selection have been performed for dataset"
)


# ## Check example output file to confirm that the process worked

# In[4]:


# Check output file
test_df = pd.read_parquet(output_feature_select_file)

print(test_df.shape)
test_df.head(5)

