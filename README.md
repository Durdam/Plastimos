# C2PLIA
## Cancer Plasticity Profiling by Live-cell Image Analysis

![image_alt](https://github.com/Durdam/C2PLIA/blob/700197d24e83af85f2ec3ed05f893e6c74cb2ba1/Images/Pipeline_framework_image.png)

We developed a semi-automated imaging analysis pipeline: Cancer Plasticity Profiling by Live-cell Image Analysis (C2PLIA). This pipeline utilizes time-series images to track live cells through deep learning-based segmentation and employs a greedy algorithm to map cell trajectories, enabling the extraction of cellular phenotypic features. These features are used to study the EMT-mediated plasticity state of cells in response to the well-known EMT-inducing factors TGF-β and EGF and TGF-β. The pipeline assigns a "plasticity index" based on various parameters, including motility, morphology, and proliferation, thus providing a quantitative measure for estimating plasticity.

## Below is a simplified workflow of the steps in this pipeline:

![image_alt](https://github.com/Durdam/C2PLIA/blob/2459d9d60497a2d9eb586a501697a549d6a026c3/Images/workflow_image.png)

Steps to implement this pipeline in R:
## Using the Code
### Requirements
This code has been developed under R 4.4.1, and Ubuntu 20.04.6 LTS.
### Installation
In addition to the above libraries, the following R packages must be installed to run the pipeline:
```shell
git clone https://github.com/Durdam/C2PLIA.git
cd C2PLIA
Rscript 01_library_installation.R
```
### Cell tracking pipeline input
To run the tracking pipeline, the pipeline requres segmented greyscale masks of each frame. The current implementation can be replicated easily by running the scripts in "Run_Tracking_Pipeline" directory.
To use the scripts for your own cell lines, replicate the scripts present for MCF7 or MDA_MB231 cell lines in "Run_Tracking_Pipeline"

### Generate segmented masks for live cell images using Cellpose
Firstly, install the cellpose and create a virtual environment. Then use the cellpose pretrained models to segment the images and generate masks to be used as input to the tracking pipeline.
Information to install cellpose and create a conda environment can be found here: (https://cellpose.readthedocs.io/en/latest/installation.html)
```shell
python -m cellpose --dir ~/path/to/images/directory/ --pretrained_model cyto3 --diameter 45 --flow_threshold 0.5 --cellprob_threshold 0 --min_size 15 --no_npy --save_tif --verbose
```
Generate masks for each condition/cell line in separate folder as done in the current implemented pipeline. See example of MCF7 directory in "Run_Tracking_Pipeline"



