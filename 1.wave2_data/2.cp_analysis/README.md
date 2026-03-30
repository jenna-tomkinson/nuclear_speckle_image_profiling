# CellProfiler segmentation and feature extraction

In this module, we perform single nuclei segmentation using the DAPI channel.
Using the nuclei masks generated, features are extracted per nuclei across each channel.
Features extracted include:

- Intensity
- Granularity
- Texture
- Correlation
- Radial distribution
- Area and shape

## Perform analysis

To run the CellProfiler analysis pipeline, run the command below:

```bash
# Make sure that you are in the 2.cp_analysis folder in terminal
source analysis.sh
```
 