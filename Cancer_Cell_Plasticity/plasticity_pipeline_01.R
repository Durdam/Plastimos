
# Description: Calculate and plot features required to compute the Plasticity_Index equation

library("ggplot2")
library("xlsx")
library("here")

# Set working directory to the root of the project
workdir = here()
setwd(workdir)

# make folder structures
if(!file.exists(paste0(workdir,"RData"))) dir.create(paste0(workdir,"RData"))
if(!file.exists(paste0(workdir,"output"))) dir.create(paste0(workdir,"output"))
if(!file.exists(paste0(workdir,"images"))) dir.create(paste0(workdir,"images"))

# set file path according to the output of the Tracking Pipeline or copy all required files in "./Cancer_Cell_Plasticity/RData/"
load("./RData/MCF7_All_Dataplot.RData")
load("./RData/MB231_All_Dataplot.RData")
load("./RData/MCF7_px_distance_all.RData")
load("./RData/MB231_px_distance_all.RData")

nrep_MCF7_ctrl = length(dataplot_MCF7_ctrl)
nrep_MCF7_EGF = length(dataplot_MCF7_EGF)
nrep_MCF7_TGF = length(dataplot_MCF7_TGF)

nrep_MB231_ctrl = length(dataplot_MDA_MB231_ctrl)
nrep_MB231_EGF = length(dataplot_MDA_MB231_EGF)
nrep_MB231_TGF = length(dataplot_MDA_MB231_TGF)

##################################################################
############################ MCF7_ctrl ###########################
##################################################################
df_MCF7_ctrl = pxdist_MCF7_ctrl= NULL
for (i in 1:nrep_MCF7_ctrl) {
  
  nobj = length(dataplot_MCF7_ctrl[[i]])
  
  pxdist_MCF7_ctrl = rbind(pxdist_MCF7_ctrl, px_distance_MCF7_ctrl[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MCF7_ctrl[[i]][[j]]) < 7) next
    df = dataplot_MCF7_ctrl[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MCF7_ctrl = rbind(df_MCF7_ctrl, df)
  }
}

pxdist_MCF7_ctrl = na.omit(pxdist_MCF7_ctrl)
pxdist_MCF7_ctrl = pxdist_MCF7_ctrl[order(pxdist_MCF7_ctrl$speed, decreasing = T), ]
pxdist_MCF7_ctrl$RepObj = paste0(pxdist_MCF7_ctrl$replicate,"_", pxdist_MCF7_ctrl$object)

df_MCF7_ctrl$AvgSpeed = df_MCF7_ctrl$TotalDistance = NA
df_MCF7_ctrl$RepObj = paste0(df_MCF7_ctrl$nReplicate, "_", df_MCF7_ctrl$nObject)

ix1 = which(df_MCF7_ctrl$RepObj %in% pxdist_MCF7_ctrl$RepObj)
ix2 = match(df_MCF7_ctrl$RepObj, pxdist_MCF7_ctrl$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MCF7_ctrl$RepObj[ix1] == pxdist_MCF7_ctrl$RepObj[ix3])) print("Perfect match in same order")

df_MCF7_ctrl$AvgSpeed[ix1] = pxdist_MCF7_ctrl$speed[ix3]
df_MCF7_ctrl$TotalDistance[ix1] = pxdist_MCF7_ctrl$distance[ix3]
df_MCF7_ctrl$Displacement[ix1] = pxdist_MCF7_ctrl$displacement[ix3]
df_MCF7_ctrl = df_MCF7_ctrl[order(df_MCF7_ctrl$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MCF7_ctrl$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MCF7_ctrl$Top20q_byAvgSpeed = "no"
df_MCF7_ctrl$Top20q_byAvgSpeed[which(df_MCF7_ctrl$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MCF7_ctrl, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_ctrl") + theme_bw()
p
ggsave(filename = "./images/MCF7_ctrl_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MCF7_ctrl[which(df_MCF7_ctrl$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_ctrl top 20%") + theme_bw()
p
ggsave(filename = "./images/MCF7_ctrl_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_ctrl: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_ctrl_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_ctrl: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_ctrl_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

##################################################################
############################ MCF7_EGF ############################
##################################################################
df_MCF7_EGF = pxdist_MCF7_EGF= NULL
for (i in 1:nrep_MCF7_EGF) {
  
  nobj = length(dataplot_MCF7_EGF[[i]])
  
  pxdist_MCF7_EGF = rbind(pxdist_MCF7_EGF, px_distance_MCF7_EGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MCF7_EGF[[i]][[j]]) < 7) next
    df = dataplot_MCF7_EGF[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MCF7_EGF = rbind(df_MCF7_EGF, df)
  }
}

pxdist_MCF7_EGF = na.omit(pxdist_MCF7_EGF)
pxdist_MCF7_EGF = pxdist_MCF7_EGF[order(pxdist_MCF7_EGF$speed, decreasing = T), ]
pxdist_MCF7_EGF$RepObj = paste0(pxdist_MCF7_EGF$replicate,"_", pxdist_MCF7_EGF$object)

df_MCF7_EGF$AvgSpeed = df_MCF7_EGF$TotalDistance = NA
df_MCF7_EGF$RepObj = paste0(df_MCF7_EGF$nReplicate, "_", df_MCF7_EGF$nObject)

ix1 = which(df_MCF7_EGF$RepObj %in% pxdist_MCF7_EGF$RepObj)
ix2 = match(df_MCF7_EGF$RepObj, pxdist_MCF7_EGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MCF7_EGF$RepObj[ix1] == pxdist_MCF7_EGF$RepObj[ix3])) print("Perfect match in same order")

df_MCF7_EGF$AvgSpeed[ix1] = pxdist_MCF7_EGF$speed[ix3]
df_MCF7_EGF$TotalDistance[ix1] = pxdist_MCF7_EGF$distance[ix3]
df_MCF7_EGF$Displacement[ix1] = pxdist_MCF7_EGF$displacement[ix3]
df_MCF7_EGF = df_MCF7_EGF[order(df_MCF7_EGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MCF7_EGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MCF7_EGF$Top20q_byAvgSpeed = "no"
df_MCF7_EGF$Top20q_byAvgSpeed[which(df_MCF7_EGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MCF7_EGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_EGF") + theme_bw()
p
ggsave(filename = "./images/MCF7_EGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MCF7_EGF[which(df_MCF7_EGF$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_EGF top 20%") + theme_bw()
p
ggsave(filename = "./images/MCF7_EGF_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_EGF: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_EGF_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_EGF: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_EGF_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

############################################################
###################### MCF7_TGF ############################
############################################################
df_MCF7_TGF = pxdist_MCF7_TGF= NULL
for (i in 1:nrep_MCF7_TGF) {
  
  nobj = length(dataplot_MCF7_TGF[[i]])
  
  pxdist_MCF7_TGF = rbind(pxdist_MCF7_TGF, px_distance_MCF7_TGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MCF7_TGF[[i]][[j]]) < 7) next
    df = dataplot_MCF7_TGF[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MCF7_TGF = rbind(df_MCF7_TGF, df)
  }
}

pxdist_MCF7_TGF = na.omit(pxdist_MCF7_TGF)
pxdist_MCF7_TGF = pxdist_MCF7_TGF[order(pxdist_MCF7_TGF$speed, decreasing = T), ]
pxdist_MCF7_TGF$RepObj = paste0(pxdist_MCF7_TGF$replicate,"_", pxdist_MCF7_TGF$object)

df_MCF7_TGF$AvgSpeed = df_MCF7_TGF$TotalDistance = NA
df_MCF7_TGF$RepObj = paste0(df_MCF7_TGF$nReplicate, "_", df_MCF7_TGF$nObject)

ix1 = which(df_MCF7_TGF$RepObj %in% pxdist_MCF7_TGF$RepObj)
ix2 = match(df_MCF7_TGF$RepObj, pxdist_MCF7_TGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MCF7_TGF$RepObj[ix1] == pxdist_MCF7_TGF$RepObj[ix3])) print("Perfect match in same order")

df_MCF7_TGF$AvgSpeed[ix1] = pxdist_MCF7_TGF$speed[ix3]
df_MCF7_TGF$TotalDistance[ix1] = pxdist_MCF7_TGF$distance[ix3]
df_MCF7_TGF$Displacement[ix1] = pxdist_MCF7_TGF$displacement[ix3]
df_MCF7_TGF = df_MCF7_TGF[order(df_MCF7_TGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MCF7_TGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MCF7_TGF$Top20q_byAvgSpeed = "no"
df_MCF7_TGF$Top20q_byAvgSpeed[which(df_MCF7_TGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MCF7_TGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_TGF") + theme_bw()
p
ggsave(filename = "./images/MCF7_TGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MCF7_TGF[which(df_MCF7_TGF$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MCF7_TGF top 20%") + theme_bw()
p
ggsave(filename = "./images/MCF7_TGF_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_TGF: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_TGF_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MCF7_TGF: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MCF7_TGF_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

##################################################################
######################### MB231_ctrl #############################
##################################################################
df_MB231_ctrl = pxdist_MB231_ctrl= NULL
for (i in 1:nrep_MB231_ctrl) {
  
  nobj = length(dataplot_MDA_MB231_ctrl[[i]])
  pxdist_MB231_ctrl = rbind(pxdist_MB231_ctrl, px_distance_MDA_MB231_ctrl[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MDA_MB231_ctrl[[i]][[j]]) < 7) next
    df = dataplot_MDA_MB231_ctrl[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MB231_ctrl = rbind(df_MB231_ctrl, df)
  }
}

pxdist_MB231_ctrl = na.omit(pxdist_MB231_ctrl)
pxdist_MB231_ctrl = pxdist_MB231_ctrl[order(pxdist_MB231_ctrl$speed, decreasing = T), ]
pxdist_MB231_ctrl$RepObj = paste0(pxdist_MB231_ctrl$replicate,"_", pxdist_MB231_ctrl$object)

df_MB231_ctrl$AvgSpeed = df_MB231_ctrl$TotalDistance = NA
df_MB231_ctrl$RepObj = paste0(df_MB231_ctrl$nReplicate, "_", df_MB231_ctrl$nObject)

ix1 = which(df_MB231_ctrl$RepObj %in% pxdist_MB231_ctrl$RepObj)
ix2 = match(df_MB231_ctrl$RepObj, pxdist_MB231_ctrl$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MB231_ctrl$RepObj[ix1] == pxdist_MB231_ctrl$RepObj[ix3])) print("Perfect match in same order")

df_MB231_ctrl$AvgSpeed[ix1] = pxdist_MB231_ctrl$speed[ix3]
df_MB231_ctrl$TotalDistance[ix1] = pxdist_MB231_ctrl$distance[ix3]
df_MB231_ctrl$Displacement[ix1] = pxdist_MB231_ctrl$displacement[ix3]
df_MB231_ctrl = df_MB231_ctrl[order(df_MB231_ctrl$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MB231_ctrl$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MB231_ctrl$Top20q_byAvgSpeed = "no"
df_MB231_ctrl$Top20q_byAvgSpeed[which(df_MB231_ctrl$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MB231_ctrl, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_ctrl") + theme_bw()
# p
ggsave(filename = "./images/MB231_ctrl_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MB231_ctrl[which(df_MB231_ctrl$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_ctrl top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_ctrl_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_ctrl: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_ctrl_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_ctrl: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_ctrl_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

##################################################################
######################### MB231_EGF ##############################
##################################################################
df_MB231_EGF = pxdist_MB231_EGF= NULL
for (i in 1:nrep_MB231_EGF) {
  
  nobj = length(dataplot_MDA_MB231_EGF[[i]])
  pxdist_MB231_EGF = rbind(pxdist_MB231_EGF, px_distance_MDA_MB231_EGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MDA_MB231_EGF[[i]][[j]]) < 7) next
    df = dataplot_MDA_MB231_EGF[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MB231_EGF = rbind(df_MB231_EGF, df)
  }
}

pxdist_MB231_EGF = na.omit(pxdist_MB231_EGF)
pxdist_MB231_EGF = pxdist_MB231_EGF[order(pxdist_MB231_EGF$speed, decreasing = T), ]
pxdist_MB231_EGF$RepObj = paste0(pxdist_MB231_EGF$replicate,"_", pxdist_MB231_EGF$object)

df_MB231_EGF$AvgSpeed = df_MB231_EGF$TotalDistance = NA
df_MB231_EGF$RepObj = paste0(df_MB231_EGF$nReplicate, "_", df_MB231_EGF$nObject)

ix1 = which(df_MB231_EGF$RepObj %in% pxdist_MB231_EGF$RepObj)
ix2 = match(df_MB231_EGF$RepObj, pxdist_MB231_EGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MB231_EGF$RepObj[ix1] == pxdist_MB231_EGF$RepObj[ix3])) print("Perfect match in same order")

df_MB231_EGF$AvgSpeed[ix1] = pxdist_MB231_EGF$speed[ix3]
df_MB231_EGF$TotalDistance[ix1] = pxdist_MB231_EGF$distance[ix3]
df_MB231_EGF$Displacement[ix1] = pxdist_MB231_EGF$displacement[ix3]
df_MB231_EGF = df_MB231_EGF[order(df_MB231_EGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MB231_EGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MB231_EGF$Top20q_byAvgSpeed = "no"
df_MB231_EGF$Top20q_byAvgSpeed[which(df_MB231_EGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MB231_EGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_EGF") + theme_bw()
p
ggsave(filename = "./images/MB231_EGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MB231_EGF[which(df_MB231_EGF$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_EGF top 20%") + theme_bw()
p
ggsave(filename = "./images/MB231_EGF_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_EGF: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_EGF_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_EGF: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_EGF_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

##################################################################
######################### MB231_EGF ##############################
##################################################################
df_MB231_TGF = pxdist_MB231_TGF= NULL
for (i in 1:nrep_MB231_TGF) {
  
  nobj = length(dataplot_MDA_MB231_TGF[[i]])
  pxdist_MB231_TGF = rbind(pxdist_MB231_TGF, px_distance_MDA_MB231_TGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_MDA_MB231_TGF[[i]][[j]]) < 7) next
    df = dataplot_MDA_MB231_TGF[[i]][[j]]
    
    # per frame speed in micrometer/hour (Time for 1 frame = 10 minutes i.e. in hours 1/6 = 0.1666667)
    df$SpeedPerFrame = c(0,
                         (sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46) / 0.1666667)
    df$DistancePerFrame = c(0,
                            sqrt(rowSums(diff(as.matrix(df[,c("m.cx","m.cy")]))^2))*0.46)
    
    # FormFactor: (2D only) Calculated as 4*π*Area/Perimeter2. Equals 1 for a perfectly circular object.
    df$FormFactor = 4*pi*df$s.area/(df$s.perimeter)^2
    
    # Calculate pearson/spearman correlation
    df$corr_Ecc_SpeedPf_pearson = cor(df$m.eccentricity, df$SpeedPerFrame, method = "pearson")
    df$corr_Ecc_SpeedPf_spearman = cor(df$m.eccentricity, df$SpeedPerFrame, method = "spearman")
    df = cbind(df, nReplicate = i, nObject = j)
    df_MB231_TGF = rbind(df_MB231_TGF, df)
  }
}

pxdist_MB231_TGF = na.omit(pxdist_MB231_TGF)
pxdist_MB231_TGF = pxdist_MB231_TGF[order(pxdist_MB231_TGF$speed, decreasing = T), ]
pxdist_MB231_TGF$RepObj = paste0(pxdist_MB231_TGF$replicate,"_", pxdist_MB231_TGF$object)

df_MB231_TGF$AvgSpeed = df_MB231_TGF$TotalDistance = NA
df_MB231_TGF$RepObj = paste0(df_MB231_TGF$nReplicate, "_", df_MB231_TGF$nObject)

ix1 = which(df_MB231_TGF$RepObj %in% pxdist_MB231_TGF$RepObj)
ix2 = match(df_MB231_TGF$RepObj, pxdist_MB231_TGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_MB231_TGF$RepObj[ix1] == pxdist_MB231_TGF$RepObj[ix3])) print("Perfect match in same order")

df_MB231_TGF$AvgSpeed[ix1] = pxdist_MB231_TGF$speed[ix3]
df_MB231_TGF$TotalDistance[ix1] = pxdist_MB231_TGF$distance[ix3]
df_MB231_TGF$Displacement[ix1] = pxdist_MB231_TGF$displacement[ix3]
df_MB231_TGF = df_MB231_TGF[order(df_MB231_TGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_MB231_TGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_MB231_TGF$Top20q_byAvgSpeed = "no"
df_MB231_TGF$Top20q_byAvgSpeed[which(df_MB231_TGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_MB231_TGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_TGF") + theme_bw()
p
ggsave(filename = "./images/MB231_TGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
dftopq = df_MB231_TGF[which(df_MB231_TGF$AvgSpeed >= top20q), ]
p = ggplot(dftopq, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("MB231_TGF top 20%") + theme_bw()
p
ggsave(filename = "./images/MB231_TGF_top20_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=TotalDistance, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_TGF: Total distance top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_TGF_top20_dist_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

# top q
p = ggplot(dftopq, aes(x=Displacement, y=m.eccentricity, color=AvgSpeed)) +
  geom_point() + ggtitle("MB231_TGF: Displacement top 20%") + theme_bw()
# p
ggsave(filename = "./images/MB231_TGF_top20_displacement_colorSpeed.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

#########################################################################
########################## Save the results #############################
#########################################################################

df_MCF7_ctrl$Cell_Line = "MCF7_ctrl"
df_MCF7_EGF$Cell_Line = "MCF7_EGF"
df_MCF7_TGF$Cell_Line = "MCF7_TGF"
df_MB231_ctrl$Cell_Line = "MB231_ctrl"
df_MB231_EGF$Cell_Line = "MB231_EGF"
df_MB231_TGF$Cell_Line = "MB231_TGF"

dfAll = rbind(df_MCF7_ctrl, df_MCF7_EGF, df_MCF7_TGF,
              df_MB231_ctrl, df_MB231_EGF, df_MB231_TGF)

dfAll$CL_RepObj = paste0(dfAll$Cell_Line,"_",dfAll$RepObj)

########### Min-Max scaling of FormFactor #############
scale_values = function(x){(x-min(x))/(max(x)-min(x))}
dfAll$FormFactor = scale_values(dfAll$FormFactor)

save(df_MCF7_ctrl, df_MCF7_EGF, df_MCF7_TGF, df_MB231_ctrl, df_MB231_EGF, df_MB231_TGF, dfAll,
     file = "./RData/All_Processed_MCF7_MB231.RData")

write.table(dfAll, file = "./output/All_Data_PostProcessed.tsv", col.names = TRUE, row.names = FALSE, sep = "\t")

# write.xlsx(dfAll, file = "./output/All_Data_PostProcessed.xlsx",
#            sheetName = "Final_processed_output", append = FALSE, row.names = FALSE)
