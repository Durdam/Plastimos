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
Generate masks for each condition/cell line in separate folder as done in the given example directory **"ImageData"**. The first 1 hour of imaging data along with segmented masks is provided in the **"ImageData"** for different conditions.

### Cell tracking pipeline
To run the tracking pipeline, segmented greyscale masks of each frame are required. To implement this on your own imaging data run the scripts in **"Run_Tracking_Pipeline"** directory.
Recommended to run it in RStudio for better visualization and inspection.
```shell
cd Run_Tracking_Pipeline
Rscript 01_forward_labeling.R
Rscript 02_feature_calculation.R
```

### Calculate Plasticity Index for each cell
To compute the **Plasticity_Index**, visualize the morphology, motility and proliferation features for different conditions and to run the tSNE algorithm on the cell features use the scripts in **"Cancer_Cell_Plasticity"** directory. Recommended to run scripts in RStudio for better visualization and inspection.
```shell
cd Cancer_Cell_Plasticity
Rscript plasticity_pipeline_01.R
Rscript plasticity_pipeline_02.R
Rscript plasticity_pipeline_03.R
Rscript plasticity_pipeline_04.R
Rscript plasticity_pipeline_05.R
```






