
# Description: Cell tracking main script. Input data: segmented masks as greyscale.

library("scales")
library("EBImage")
library("foreach")
library("doParallel")
library("here")

outputfilename = "Cell_Line_x_All_Dataplot.RData"

effZ = (.Machine$double.eps)^0.95

# Set working directory to the root of the project
workdir = here()
setwd(workdir)

# make folder structures
if(!file.exists(paste0(workdir,"/RData"))) dir.create(paste0(workdir,"/RData"))
if(!file.exists(paste0(workdir,"/output"))) dir.create(paste0(workdir,"/output"))
if(!file.exists(paste0(workdir,"/images"))) dir.create(paste0(workdir,"/images"))

# source package R functions from the primary cell tracking functions
package.path = "../Cell_Tracking_Primary_Functions"
Rfiles = list.files(path=package.path, pattern="[[:print:]]*.R",full.names=TRUE)
for(fR in Rfiles) source(fR)

## The cell line with different treatments = replicates with imaging data
Cell_Line_x_ctrl = c("A1","A2")
Cell_Line_x_EGF = c("A1","A2")
Cell_Line_x_TGF = c("A1","A2")

## Path to folder containing imaging data
file_Cell_Line_x_ctrl = paste0("../ImageData/ctrl/", Cell_Line_x_ctrl,"/")
file_Cell_Line_x_EGF = paste0("../ImageData/egf/", Cell_Line_x_EGF,"/")
file_Cell_Line_x_TGF = paste0("../ImageData/tgf/", Cell_Line_x_TGF,"/")

# Number of frames to track the cells: 150 frames = 25 hrs
nframekeep = 150

######################################### USER INPUT ENDS ########################################################
cat("Starting cell tracking workflow....\n")

#############################
##### Cell_Line_x_ctrl ######
#############################
dataplot_Cell_Line_x_ctrl = list()
for (j in 1:length(file_Cell_Line_x_ctrl)) {
  
  iname = list.files(file_Cell_Line_x_ctrl[j], pattern = "cp_masks.tif")
  if(length(iname) > nframekeep) iname = tail(iname, nframekeep)
  iname = data.frame(iname)
  colnames(iname) = "imgname"
  nimg = nrow(iname)
  cellmask = feat = list()
  
  for(i in 1:nimg) {
    
    imgpath = paste0(file_Cell_Line_x_ctrl[j], iname[i,])
    tifmask = readImage(imgpath)
    unimask = sort(unique(c(tifmask)))
    off = ifelse(unimask[1] == 0,-1,0)
    intmask = tifmask
    intmask[] = 0
    
    for (k in 1:length(unimask)) {
      intmask[tifmask == unimask[k]] = k+off
    }
    ################
    testmask = splitMultiMask(intmask) ## mask for every single cell in a particular frame
    cellmask[[i]] = testmask
    ################
    
    #### xy cordinate calculation
    feat_moment = computeFeatures.moment(intmask)
    feat_shape = computeFeatures.shape(intmask)
    feat[[i]] = cbind(Index = as.numeric(rownames(feat_moment)), feat_moment, feat_shape)
  } 
  
  celllabels = Labeling(cellmask, multimask = FALSE, ncoressys = 0)

  ## matching the labels with cell coordinates
  dataxy = list()
  datacp = list()
  ix1 = ix2 = NA
  for (i in 1:nimg) {
    celllabels[[i]] = cbind(celllabels[[i]], frame = rep(i, nrow(celllabels[[i]])))
    dataxy[[i]] = cbind(celllabels[[i]], feat[[i]]) 
  }
  
  datacp = Reduce(rbind, dataxy) # convert list to matrix
  datacp = as.data.frame(datacp)
  
  # list of tracked cells for all frames
  datacpList = split(datacp, f = datacp$label)
  # ngroups = length(datacpList)
  dataplot_Cell_Line_x_ctrl[[j]] = datacpList
  
  ## status
  print(paste0("Current running: ",file_Cell_Line_x_ctrl[j]))
}

############################
##### Cell_Line_x_EGF ######
############################
dataplot_Cell_Line_x_EGF = list()
for (j in 1:length(file_Cell_Line_x_EGF)) {
  
  iname = list.files(file_Cell_Line_x_EGF[j], pattern = "cp_masks.tif")
  if(length(iname) > nframekeep) iname = tail(iname, nframekeep)
  iname = data.frame(iname)
  colnames(iname) = "imgname"
  nimg = nrow(iname)
  cellmask = feat = list() 
  
  for(i in 1:nimg) {
    
    imgpath = paste0(file_Cell_Line_x_EGF[j], iname[i,])
    tifmask = readImage(imgpath)
    unimask = sort(unique(c(tifmask)))
    off = ifelse(unimask[1] == 0,-1,0)
    intmask = tifmask
    intmask[] = 0
    
    for (k in 1:length(unimask)) {
      intmask[tifmask == unimask[k]] = k+off
    }
    ################
    testmask = splitMultiMask(intmask) ## mask for every single cell in a particular frame
    cellmask[[i]] = testmask
    ################
    
    #### xy cordinate calculation
    feat_moment = computeFeatures.moment(intmask)
    feat_shape = computeFeatures.shape(intmask)
    feat[[i]] = cbind(Index = as.numeric(rownames(feat_moment)), feat_moment, feat_shape)
  } 
  
  celllabels = Labeling(cellmask, multimask = FALSE, ncoressys = 0)
  
  ## matching the labels with cell coordinates
  dataxy = list()
  datacp = list()
  ix1 = ix2 = NA
  for (i in 1:nimg) {
    celllabels[[i]] = cbind(celllabels[[i]], frame = rep(i, nrow(celllabels[[i]])))
    dataxy[[i]] = cbind(celllabels[[i]], feat[[i]]) 
  }
  
  datacp = Reduce(rbind, dataxy) # convert list to matrix
  datacp = as.data.frame(datacp)
  
  # list of tracked cells for all frames
  datacpList = split(datacp, f = datacp$label)
  # ngroups = length(datacpList)
  dataplot_Cell_Line_x_EGF[[j]] = datacpList
  
  ## status
  print(paste0("Current running: ",file_Cell_Line_x_EGF[j]))
}

############################
##### Cell_Line_x_TGF ######
############################
dataplot_Cell_Line_x_TGF = list()
for (j in 1:length(file_Cell_Line_x_TGF)) {
  
  iname = list.files(file_Cell_Line_x_TGF[j], pattern = "cp_masks.tif")
  if(length(iname) > nframekeep) iname = tail(iname, nframekeep)
  iname = data.frame(iname)
  colnames(iname) = "imgname"
  nimg = nrow(iname)
  cellmask = feat = list() 
  
  for(i in 1:nimg) {
    
    imgpath = paste0(file_Cell_Line_x_TGF[j], iname[i,])
    tifmask = readImage(imgpath)
    unimask = sort(unique(c(tifmask)))
    off = ifelse(unimask[1] == 0,-1,0)
    intmask = tifmask
    intmask[] = 0
    
    for (k in 1:length(unimask)) {
      intmask[tifmask == unimask[k]] = k+off
    }
    ################
    testmask = splitMultiMask(intmask) ## mask for every single cell in a particular frame
    cellmask[[i]] = testmask
    ################
    
    #### xy cordinate calculation
    feat_moment = computeFeatures.moment(intmask)
    feat_shape = computeFeatures.shape(intmask)
    feat[[i]] = cbind(Index = as.numeric(rownames(feat_moment)), feat_moment, feat_shape)
  } 
  
  celllabels = Labeling(cellmask, multimask = FALSE, ncoressys = 0)
  
  ## matching the labels with cell coordinates
  dataxy = list()
  datacp = list()
  ix1 = ix2 = NA
  for (i in 1:nimg) {
    celllabels[[i]] = cbind(celllabels[[i]], frame = rep(i, nrow(celllabels[[i]])))
    dataxy[[i]] = cbind(celllabels[[i]], feat[[i]]) 
  }
  
  datacp = Reduce(rbind, dataxy) # convert list to matrix
  datacp = as.data.frame(datacp)
  
  # list of tracked cells for all frames
  datacpList = split(datacp, f = datacp$label)
  # ngroups = length(datacpList)
  dataplot_Cell_Line_x_TGF[[j]] = datacpList
  
  ## status
  print(paste0("Current running: ",file_Cell_Line_x_TGF[j]))
}


########################## Save Cell Tracking Results ###################################
save(dataplot_Cell_Line_x_ctrl, dataplot_Cell_Line_x_EGF, dataplot_Cell_Line_x_TGF,
     file = paste0("./RData/", outputfilename))
#########################################################################################

if(file.exists(paste0(workdir,"/RData/",outputfilename))) {
cat(paste0("✅ Script completed successfully!\nResults saved to './RData/",outputfilename,"'.\n"))
} else cat("❌ Script failed. Please check if input files are provided correctly.\n")

