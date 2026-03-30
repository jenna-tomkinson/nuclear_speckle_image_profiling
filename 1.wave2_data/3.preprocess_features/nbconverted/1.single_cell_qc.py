#!/usr/bin/env python
# coding: utf-8

# # Perform single-cell quality control
# 
# In this notebook, we perform single-cell quality control using coSMicQC.
# To filter the single-cells, the default method is z-score to find outliers using the values from only one feature at a time. 
# We use features from the AreaShape and Intensity modules to assess the quality of the segmented single-cells:
# 
# ### Assessing poor nuclei segmentation
# 
# We use the following features either in combination or alone when finding technical segmentation errors.
# To be specific, we use:
# 
# | Measurement(s)           | Target                                                           | Notes                                                                                                                                                        |
# |--------------------------|------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
# | Nuclei Soldity + Mass Displacement in DAPI channel | Mis-segmented nuclei from debris, artifacts or over-segmented clusters | We use these two features in combination where one finds the indented/non-smooth nuclei segmentations and the finds segmentations where the centeroid based on intensity is highly shifted (meaning multiple nuclei or segmented artifacts).                                           |
# | Nuclei Area      | Under-segmented nuclei   | We use this by itself to find abnormally small nuclei which related to under-segmentation from the algorithm.                                                  |
# 

# ## Import libraries

# In[1]:


import pathlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from cytodataframe import CytoDataFrame
from cosmicqc import find_outliers


# ## Set paths and variables

# In[2]:


# Directory with data
data_dir = pathlib.Path("./data/converted_profiles/")

# Directory to save cleaned data
qc_dir = pathlib.Path("./data/single_cell_qc/")
qc_dir.mkdir(exist_ok=True)

# Directory to save qc figures
qc_fig_dir = pathlib.Path("./qc_figures")
qc_fig_dir.mkdir(exist_ok=True)

# Create an empty dictionary to store data frames for each plate
all_qc_data_frames = {}

# metadata columns to include in output data frame
metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Position",
    "Image_Metadata_Best_Z",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
    "Nuclei_AreaShape_BoundingBoxMaximum_X",
    "Nuclei_AreaShape_BoundingBoxMaximum_Y",
    "Nuclei_AreaShape_BoundingBoxMinimum_X",
    "Nuclei_AreaShape_BoundingBoxMinimum_Y",
    "Image_PathName_DAPI",
    "Image_FileName_DAPI",
]


# ## Concat all plates together to assess quality control

# In[3]:


# Only process the files that are in the plate names list
plate_path = pathlib.Path(f"{data_dir}/u2os_per_nuclei.parquet")

# Read the parquet file into a DataFrame
qc_df = pd.read_parquet(plate_path)

print(qc_df.shape)
qc_df.head()


# ## Evaluate the distributions of the features across plates to determine best features to use for QC

# In[4]:


# Set the style
sns.set_style("whitegrid")

# Create subplots with 2 rows and 2 columns
fig, axes = plt.subplots(1, 3, figsize=(18, 6))  # first is width and second is height

# Flatten the axes array for easy indexing
axes = axes.flatten()

# Create kdeplot for Nuclei_AreaShape_Area
sns.kdeplot(
    data=qc_df,
    x="Nuclei_AreaShape_Area",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[0],  # Set the first subplot
)
axes[0].set_title("Area")

# Create kdeplot for FormFactor
sns.kdeplot(
    data=qc_df,
    x="Nuclei_AreaShape_FormFactor",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[1],
    legend=False,
)
axes[1].set_title("Form Factor")

# Create kdeplot for MassDisplacement
sns.kdeplot(
    data=qc_df,
    x="Nuclei_Intensity_MassDisplacement_DAPI",
    hue="Image_Metadata_Plate",
    palette="viridis",
    fill=True,
    common_norm=False,
    ax=axes[2],
    legend=False,
)
axes[2].set_title("Mass Displacement")

plt.suptitle(
    "Distribution of each quality control feature across all plates",
    fontsize=16,
)
plt.tight_layout(rect=[0, 0, 1.0, 1.0])  # Adjust subplot layout

# Save figure
plt.savefig(pathlib.Path(f"{qc_fig_dir}/QC_features_dist_plot.png"), dpi=500)

plt.show()


# ## Use Nuclei Solidity and Mass Displacement to detect improperly segmented nuclei (debris, artifacts, over-segmentation)

# In[5]:


outline_to_orig_mapping = {
    rf"{record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}_{record['Image_Metadata_Position']}_ch01_{record['Image_Metadata_Best_Z']}_illumcorrect_MaskNuclei.tiff": rf"{record['Image_Metadata_Plate']}_{record['Image_Metadata_Well']}_{record['Image_Metadata_Position']}_ch\d+_{record['Image_Metadata_Best_Z']}_illumcorrect.tiff"
    for record in qc_df[
        [
            "Image_Metadata_Plate",
            "Image_Metadata_Well",
            "Image_Metadata_Position",
            "Image_Metadata_Best_Z",
        ]
    ].to_dict(orient="records")
}

next(iter(outline_to_orig_mapping.items()))


# In[6]:


# find abnormally large nuclei based on area
feature_thresholds = {
    "Nuclei_AreaShape_Solidity": -2,
    "Nuclei_Intensity_MassDisplacement_DAPI": 1,
}

missegmented_nuclei_outliers = find_outliers(
    df=qc_df, metadata_columns=metadata_columns, feature_thresholds=feature_thresholds
)

# Convert to CytoDataFrame for outline viewing
missegmented_nuclei_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(missegmented_nuclei_outliers),
    data_context_dir="../1.illumination_correction/IC_corrected_images",
    data_mask_context_dir="../2.cp_analysis/analysis_output/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
    },
)[
    [
        "Nuclei_AreaShape_Solidity",
        "Nuclei_Intensity_MassDisplacement_DAPI",
        "Image_FileName_DAPI",
    ]
]

missegmented_nuclei_outliers_cdf.sort_values(
    by="Nuclei_AreaShape_Solidity", ascending=True
).head().T


# In[7]:


# find abnormally small nuclei based on area
feature_thresholds = {
    "Nuclei_AreaShape_Area": -2,
}

small_nuclei_outliers = find_outliers(
    df=qc_df, metadata_columns=metadata_columns, feature_thresholds=feature_thresholds
)

# Convert to CytoDataFrame for outline viewing
small_nuclei_outliers_cdf = CytoDataFrame(
    data=pd.DataFrame(small_nuclei_outliers),
    data_context_dir="../1.illumination_correction/IC_corrected_images",
    data_mask_context_dir="../2.cp_analysis/analysis_output/",
    segmentation_file_regex=outline_to_orig_mapping,
    display_options={
        "center_dot": True,
        "brightness": 5,
        "outline_color": (180, 30, 180),
    },
)[
    [
        "Nuclei_AreaShape_Area",
        "Image_FileName_DAPI",
    ]
]

small_nuclei_outliers_cdf.sort_values(
    by="Nuclei_AreaShape_Area", ascending=True
).head().T


# ## Output file with metadata and columns indicating if the cell failed QC or not

# In[8]:


# Set metadata columns to include in output data frame
output_metadata_columns = [
    "Image_Metadata_Plate",
    "Image_Metadata_Well",
    "Image_Metadata_Position",
    "Image_Metadata_Best_Z",
    "Nuclei_Location_Center_X",
    "Nuclei_Location_Center_Y",
]

# Make a copy of qc_df to avoid modifying the original DataFrame
qc_df_copy = qc_df.copy()

# Add columns to indicate QC failure for missegmented and small nuclei
qc_df_copy["Metadata_cqc_failed_missegmented_nuclei"] = qc_df_copy.index.isin(
    missegmented_nuclei_outliers.index
)
qc_df_copy["Metadata_cqc_failed_small_nuclei"] = qc_df_copy.index.isin(
    small_nuclei_outliers.index
)

# Select only metadata columns and the new QC columns
qc_output_df = qc_df_copy[
    output_metadata_columns
    + [
        "Metadata_cqc_failed_missegmented_nuclei",
        "Metadata_cqc_failed_small_nuclei",
    ]
]

# Save the resulting DataFrame to a parquet file
qc_output_df.to_parquet(qc_dir / "qc.parquet", index=False)

# Show the first few rows of the output DataFrame
print(
    "Number of nuclei that failed QC in either condition:",
    qc_output_df[
        ["Metadata_cqc_failed_missegmented_nuclei", "Metadata_cqc_failed_small_nuclei"]
    ]
    .any(axis=1)
    .sum(),
)
print("Shape of the output DataFrame:", qc_output_df.shape)
qc_output_df.head()

