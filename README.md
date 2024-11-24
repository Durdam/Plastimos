# C2PLIA
## Cancer Plasticity Profiling by Live-cell Image Analysis

![image_alt](https://github.com/Durdam/C2PLIA/blob/700197d24e83af85f2ec3ed05f893e6c74cb2ba1/Images/Pipeline_framework_image.png)

We developed a semi-automated imaging analysis pipeline: **C**ancer **P**lasticity **P**rofiling by **L**ive-cell **I**mage **A**nalysis (C2PLIA). This pipeline utilizes time-series images to track live cells through deep learning-based segmentation and employs a greedy algorithm to map cell trajectories, enabling the extraction of cellular phenotypic features. These features are used to study the EMT-mediated plasticity state of cells in response to the well-known EMT-inducing factors EGF and TGF-β. The pipeline assigns a **"Plasticity Index"** based on various parameters, including motility, morphology, and proliferation, thus providing a quantitative measure for estimating plasticity. This scalable, live cell imaging based framework offers a powerful tool for quantifying EMT-mediated plasticity and can be applied to high-throughput drug screening.

## Below is a simplified workflow of the sequence of steps in this pipeline:

![image_alt](https://github.com/Durdam/C2PLIA/blob/2459d9d60497a2d9eb586a501697a549d6a026c3/Images/workflow_image.png)

## Steps to implement this pipeline in R:
### Requirements
This code has been developed under **R 4.4.1**, and **Ubuntu 20.04.6**.
### Installation
In addition to the above libraries, the following R packages must be installed to run the pipeline:
```shell
git clone https://github.com/Durdam/C2PLIA.git
cd C2PLIA
Rscript 01_library_installation.R
```

### Generate segmented masks for live cell images using Cellpose
Below is an example on how you can generate segmented masks from your raw microscopy images to give input to this pipeline. 
Firstly, install the cellpose and create a virtual environment. Then use the cellpose pretrained models to segment the images and generate masks to be used as input to the tracking pipeline.
Information to install cellpose and create a conda environment can be found here: (https://cellpose.readthedocs.io/en/latest/installation.html)
**Generate masks:**
```shell
python -m cellpose --dir ~/path/to/images/directory/ --pretrained_model cyto3 --diameter 45 --flow_threshold 0.5 --cellprob_threshold 0 --min_size 15 --no_npy --save_tif --verbose
```
Generate masks for each condition/cell_line in separate folders as shown in the given example directory **"ImageData"**. The first 3 hour of imaging data along with segmented masks are provided in the **"ImageData"** for different conditions from MCF7 cell line as example input data for the pipeline.

### Cell tracking pipeline
To run the tracking pipeline, segmented greyscale masks of each frame are required. To implement this on your own imaging data run the scripts in **"Run_Tracking_Pipeline"** directory. The primary tracking functions are available in **"Cell_Tracking_Primary_Functions"** directory which are called from `01_forward_labeling.R`
```shell
cd Run_Cell_Tracking
Rscript 01_forward_labeling.R
Rscript 02_feature_calculation.R
cd ..
```

### Calculate Plasticity_Index for cells
To compute the **Plasticity_Index** along with visualizing the morphology, motility and proliferation features for different conditions and to run the tSNE algorithm on the cell features use the scripts in **"Cancer_Cell_Plasticity"** directory. Below is implementation in shell.
```shell
cd Cancer_Cell_Plasticity

# Calculate features required to compute the Plasticity_Index equation
Rscript plasticity_pipeline_01.R

# Plotting and visualizations of the extracted cell features
Rscript plasticity_pipeline_02.R

# Compute growth rate: Required parameters:
# nframe = number of frames for which you have imaging data
# nframe_gap = Time gap (in minutes) between each frame
# time_0_frame = "y" or "n". Select "y" if your imaging data has image at time = 0
Rscript plasticity_pipeline_03.R nframe=19 nframe_gap=10 time_0_frame=y

# Compute Plasticity_Index per cell
Rscript plasticity_pipeline_04.R

# Run tSNE on Plasticity_Index parameters and plotting
Rscript plasticity_pipeline_05.R
```






