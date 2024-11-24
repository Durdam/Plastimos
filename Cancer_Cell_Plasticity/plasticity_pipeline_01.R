
# Description: Calculate and plot features required to compute the Plasticity_Index equation

library("ggplot2")
library("xlsx")
library("here")

outputfilename_RData = "All_Processed_Cell_Line_x.RData"
outputfilename_output = "All_Data_PostProcessed.tsv"

# Set working directory to the root of the project
workdir = here()
setwd(workdir)

# make folder structures
if(!file.exists(paste0(workdir,"/RData"))) dir.create(paste0(workdir,"/RData"))
if(!file.exists(paste0(workdir,"/output"))) dir.create(paste0(workdir,"/output"))
if(!file.exists(paste0(workdir,"/images"))) dir.create(paste0(workdir,"/images"))

# set file path according to the output of the Tracking Pipeline or copy all required files in "./Cancer_Cell_Plasticity/RData/"
load("../Run_Cell_Tracking/RData/Cell_Line_x_All_Dataplot.RData")
load("../Run_Cell_Tracking/RData/Cell_Line_x_px_distance_all.RData")

nrep_Cell_Line_x_ctrl = length(dataplot_Cell_Line_x_ctrl)
nrep_Cell_Line_x_EGF = length(dataplot_Cell_Line_x_EGF)
nrep_Cell_Line_x_TGF = length(dataplot_Cell_Line_x_TGF)

print("Computing Plasticity Index parameters.....")

##################################################################
############################ Cell_Line_x_ctrl ####################
##################################################################
df_Cell_Line_x_ctrl = pxdist_Cell_Line_x_ctrl= NULL
for (i in 1:nrep_Cell_Line_x_ctrl) {
  
  nobj = length(dataplot_Cell_Line_x_ctrl[[i]])
  
  pxdist_Cell_Line_x_ctrl = rbind(pxdist_Cell_Line_x_ctrl, px_distance_Cell_Line_x_ctrl[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_Cell_Line_x_ctrl[[i]][[j]]) < 7) next
    df = dataplot_Cell_Line_x_ctrl[[i]][[j]]
    
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
    df_Cell_Line_x_ctrl = rbind(df_Cell_Line_x_ctrl, df)
  }
}

pxdist_Cell_Line_x_ctrl = na.omit(pxdist_Cell_Line_x_ctrl)
pxdist_Cell_Line_x_ctrl = pxdist_Cell_Line_x_ctrl[order(pxdist_Cell_Line_x_ctrl$speed, decreasing = T), ]
pxdist_Cell_Line_x_ctrl$RepObj = paste0(pxdist_Cell_Line_x_ctrl$replicate,"_", pxdist_Cell_Line_x_ctrl$object)

df_Cell_Line_x_ctrl$AvgSpeed = df_Cell_Line_x_ctrl$TotalDistance = NA
df_Cell_Line_x_ctrl$RepObj = paste0(df_Cell_Line_x_ctrl$nReplicate, "_", df_Cell_Line_x_ctrl$nObject)

ix1 = which(df_Cell_Line_x_ctrl$RepObj %in% pxdist_Cell_Line_x_ctrl$RepObj)
ix2 = match(df_Cell_Line_x_ctrl$RepObj, pxdist_Cell_Line_x_ctrl$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_Cell_Line_x_ctrl$RepObj[ix1] == pxdist_Cell_Line_x_ctrl$RepObj[ix3])) print("Perfect match in same order")

df_Cell_Line_x_ctrl$AvgSpeed[ix1] = pxdist_Cell_Line_x_ctrl$speed[ix3]
df_Cell_Line_x_ctrl$TotalDistance[ix1] = pxdist_Cell_Line_x_ctrl$distance[ix3]
df_Cell_Line_x_ctrl$Displacement[ix1] = pxdist_Cell_Line_x_ctrl$displacement[ix3]
df_Cell_Line_x_ctrl = df_Cell_Line_x_ctrl[order(df_Cell_Line_x_ctrl$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_Cell_Line_x_ctrl$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_Cell_Line_x_ctrl$Top20q_byAvgSpeed = "no"
df_Cell_Line_x_ctrl$Top20q_byAvgSpeed[which(df_Cell_Line_x_ctrl$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_Cell_Line_x_ctrl, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("Cell_Line_x_ctrl") + theme_bw()
# p
ggsave(filename = "./images/Cell_Line_x_ctrl_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

##################################################################
############################ Cell_Line_x_EGF #####################
##################################################################
df_Cell_Line_x_EGF = pxdist_Cell_Line_x_EGF= NULL
for (i in 1:nrep_Cell_Line_x_EGF) {
  
  nobj = length(dataplot_Cell_Line_x_EGF[[i]])
  
  pxdist_Cell_Line_x_EGF = rbind(pxdist_Cell_Line_x_EGF, px_distance_Cell_Line_x_EGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_Cell_Line_x_EGF[[i]][[j]]) < 7) next
    df = dataplot_Cell_Line_x_EGF[[i]][[j]]
    
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
    df_Cell_Line_x_EGF = rbind(df_Cell_Line_x_EGF, df)
  }
}

pxdist_Cell_Line_x_EGF = na.omit(pxdist_Cell_Line_x_EGF)
pxdist_Cell_Line_x_EGF = pxdist_Cell_Line_x_EGF[order(pxdist_Cell_Line_x_EGF$speed, decreasing = T), ]
pxdist_Cell_Line_x_EGF$RepObj = paste0(pxdist_Cell_Line_x_EGF$replicate,"_", pxdist_Cell_Line_x_EGF$object)

df_Cell_Line_x_EGF$AvgSpeed = df_Cell_Line_x_EGF$TotalDistance = NA
df_Cell_Line_x_EGF$RepObj = paste0(df_Cell_Line_x_EGF$nReplicate, "_", df_Cell_Line_x_EGF$nObject)

ix1 = which(df_Cell_Line_x_EGF$RepObj %in% pxdist_Cell_Line_x_EGF$RepObj)
ix2 = match(df_Cell_Line_x_EGF$RepObj, pxdist_Cell_Line_x_EGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_Cell_Line_x_EGF$RepObj[ix1] == pxdist_Cell_Line_x_EGF$RepObj[ix3])) print("Perfect match in same order")

df_Cell_Line_x_EGF$AvgSpeed[ix1] = pxdist_Cell_Line_x_EGF$speed[ix3]
df_Cell_Line_x_EGF$TotalDistance[ix1] = pxdist_Cell_Line_x_EGF$distance[ix3]
df_Cell_Line_x_EGF$Displacement[ix1] = pxdist_Cell_Line_x_EGF$displacement[ix3]
df_Cell_Line_x_EGF = df_Cell_Line_x_EGF[order(df_Cell_Line_x_EGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_Cell_Line_x_EGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_Cell_Line_x_EGF$Top20q_byAvgSpeed = "no"
df_Cell_Line_x_EGF$Top20q_byAvgSpeed[which(df_Cell_Line_x_EGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_Cell_Line_x_EGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("Cell_Line_x_EGF") + theme_bw()
# p
ggsave(filename = "./images/Cell_Line_x_EGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

############################################################
###################### Cell_Line_x_TGF #####################
############################################################
df_Cell_Line_x_TGF = pxdist_Cell_Line_x_TGF= NULL
for (i in 1:nrep_Cell_Line_x_TGF) {
  
  nobj = length(dataplot_Cell_Line_x_TGF[[i]])
  
  pxdist_Cell_Line_x_TGF = rbind(pxdist_Cell_Line_x_TGF, px_distance_Cell_Line_x_TGF[[i]])
  
  for (j in 1:nobj) {
    # Take only cells that occur in atleast 7 frames to reduce the amount of artifacts
    if(nrow(dataplot_Cell_Line_x_TGF[[i]][[j]]) < 7) next
    df = dataplot_Cell_Line_x_TGF[[i]][[j]]
    
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
    df_Cell_Line_x_TGF = rbind(df_Cell_Line_x_TGF, df)
  }
}

pxdist_Cell_Line_x_TGF = na.omit(pxdist_Cell_Line_x_TGF)
pxdist_Cell_Line_x_TGF = pxdist_Cell_Line_x_TGF[order(pxdist_Cell_Line_x_TGF$speed, decreasing = T), ]
pxdist_Cell_Line_x_TGF$RepObj = paste0(pxdist_Cell_Line_x_TGF$replicate,"_", pxdist_Cell_Line_x_TGF$object)

df_Cell_Line_x_TGF$AvgSpeed = df_Cell_Line_x_TGF$TotalDistance = NA
df_Cell_Line_x_TGF$RepObj = paste0(df_Cell_Line_x_TGF$nReplicate, "_", df_Cell_Line_x_TGF$nObject)

ix1 = which(df_Cell_Line_x_TGF$RepObj %in% pxdist_Cell_Line_x_TGF$RepObj)
ix2 = match(df_Cell_Line_x_TGF$RepObj, pxdist_Cell_Line_x_TGF$RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(df_Cell_Line_x_TGF$RepObj[ix1] == pxdist_Cell_Line_x_TGF$RepObj[ix3])) print("Perfect match in same order")

df_Cell_Line_x_TGF$AvgSpeed[ix1] = pxdist_Cell_Line_x_TGF$speed[ix3]
df_Cell_Line_x_TGF$TotalDistance[ix1] = pxdist_Cell_Line_x_TGF$distance[ix3]
df_Cell_Line_x_TGF$Displacement[ix1] = pxdist_Cell_Line_x_TGF$displacement[ix3]
df_Cell_Line_x_TGF = df_Cell_Line_x_TGF[order(df_Cell_Line_x_TGF$AvgSpeed, decreasing = TRUE), ]
top20q = quantile(df_Cell_Line_x_TGF$AvgSpeed, probs = seq(0, 1, 0.2))[5]
df_Cell_Line_x_TGF$Top20q_byAvgSpeed = "no"
df_Cell_Line_x_TGF$Top20q_byAvgSpeed[which(df_Cell_Line_x_TGF$AvgSpeed >= top20q)] = "yes"

# all cells
p = ggplot(df_Cell_Line_x_TGF, aes(x=SpeedPerFrame, y=m.eccentricity, color=SpeedPerFrame)) +
  geom_point() + ggtitle("Cell_Line_x_TGF") + theme_bw()
# p
ggsave(filename = "./images/Cell_Line_x_TGF_Speed_Ecc.png", plot = p, width = 8, height = 6, units = "in", dpi = 300)

#####################################################################################

print("Plotting Eccentricity vs SpeedperFrame parameter for all cells in different conditions.....")
print("Saving Results....")
#########################################################################
########################## Save the results #############################
#########################################################################

df_Cell_Line_x_ctrl$Cell_Line = "Cell_Line_x_ctrl"
df_Cell_Line_x_EGF$Cell_Line = "Cell_Line_x_EGF"
df_Cell_Line_x_TGF$Cell_Line = "Cell_Line_x_TGF"


dfAll = rbind(df_Cell_Line_x_ctrl, df_Cell_Line_x_EGF, df_Cell_Line_x_TGF)

dfAll$CL_RepObj = paste0(dfAll$Cell_Line,"_",dfAll$RepObj)

########### Min-Max scaling of FormFactor #############
scale_values = function(x){(x-min(x))/(max(x)-min(x))}
dfAll$FormFactor = scale_values(dfAll$FormFactor)

save(df_Cell_Line_x_ctrl, df_Cell_Line_x_EGF, df_Cell_Line_x_TGF, dfAll,
     file = paste0("./RData/",outputfilename_RData))

write.table(dfAll, file = paste0("./output/", outputfilename_output), col.names = TRUE, row.names = FALSE, sep = "\t")

# write.xlsx(dfAll, file = "./output/All_Data_PostProcessed.xlsx",
#            sheetName = "Final_processed_output", append = FALSE, row.names = FALSE)

output_files = c(paste0(workdir,"/RData/",outputfilename_RData), paste0(workdir,"/output/",outputfilename_output))
if(unique(file.exists(output_files))) {
  cat(paste0("✅ Script completed successfully!\nResults saved to './RData/",outputfilename_RData,"' and './output/",outputfilename_output,"'.\n"))
} else cat("❌ Script failed. Please check if input files are provided correctly.\n")
