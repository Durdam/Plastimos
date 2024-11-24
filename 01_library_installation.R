
# packages required to run the pipeline
# Description: Installation of the required R packages (cran + bioconductor repositories + devtools)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EBImage")

# Function to check and install packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg))
      install.packages(pkg)
    } else {
      message(paste(pkg, "is already installed."))
    }
  }
}

# List of required packages
required_packages <- c("dplyr", "ggplot2", "here", "tidyr", "scales", "colortools", "RColorBrewer", 
"plot3D", "foreach", "doParallel", "factoextra", "xlsx", "Rtsne")

# Check and install
check_and_install(required_packages)

