
# Description: Cell features calculation

outputfilename = ("Cell_Line_x_px_distance_all.RData")

cat("Current working directory is set to:", getwd(), "\n")
# Set working directory to the root of the project
workdir = getwd()
setwd(workdir)

# output of previous script: 01_forward_labelling.R
load("./RData/Cell_Line_x_All_Dataplot.RData")

#################################################################
cat("Computing morphology and motiliy features....\n")

# calculate morphology and motility features of objects/cells
######## Cell_Line_x_ctrl #########
nobject = lengths(dataplot_Cell_Line_x_ctrl)
px_distance_Cell_Line_x_ctrl = NULL

for(i in 1:length(nobject)) {
  px_distance_Cell_Line_x_ctrl[[i]] = data.frame(object = rep(NA,nobject[i]), distance = rep(NA,nobject[i]), 
                                         replicate = rep(NA,nobject[i]), speed = rep(NA,nobject[i]), 
                                         area = rep(NA,nobject[i]), perimeter = rep(NA,nobject[i]),
                                         displacement = rep(NA,nobject[i]))
  for (j in 1:nobject[i]) {
    if(nrow(dataplot_Cell_Line_x_ctrl[[i]][[j]]) < 2) next
    ix_last_frame = nrow(dataplot_Cell_Line_x_ctrl[[i]][[j]])
    px_distance_Cell_Line_x_ctrl[[i]][j,] = c(j, (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_ctrl[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46, i, 
                                      (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_ctrl[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46 / (nrow(dataplot_Cell_Line_x_ctrl[[i]][[j]]) / 6),
                                      mean(dataplot_Cell_Line_x_ctrl[[i]][[j]][,"s.area"]), mean(dataplot_Cell_Line_x_ctrl[[i]][[j]][,"s.perimeter"]),
                                      sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_ctrl[[i]][[j]] [,c("m.cx","m.cy")][c(1,ix_last_frame),]))^2))*0.46)
    
  }
}

######## Cell_Line_x_EGF #########
nobject = lengths(dataplot_Cell_Line_x_EGF)
px_distance_Cell_Line_x_EGF = NULL

for(i in 1:length(nobject)) {
  px_distance_Cell_Line_x_EGF[[i]] = data.frame(object = rep(NA,nobject[i]), distance = rep(NA,nobject[i]), 
                                         replicate = rep(NA,nobject[i]), speed = rep(NA,nobject[i]), 
                                         area = rep(NA,nobject[i]), perimeter = rep(NA,nobject[i]),
                                         displacement = rep(NA,nobject[i]))
  for (j in 1:nobject[i]) {
    if(nrow(dataplot_Cell_Line_x_EGF[[i]][[j]]) < 2) next
    ix_last_frame = nrow(dataplot_Cell_Line_x_EGF[[i]][[j]])
    px_distance_Cell_Line_x_EGF[[i]][j,] = c(j, (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_EGF[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46, i, 
                                      (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_EGF[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46 / (nrow(dataplot_Cell_Line_x_EGF[[i]][[j]]) / 6),
                                      mean(dataplot_Cell_Line_x_EGF[[i]][[j]][,"s.area"]), mean(dataplot_Cell_Line_x_EGF[[i]][[j]][,"s.perimeter"]),
                                      sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_EGF[[i]][[j]] [,c("m.cx","m.cy")][c(1,ix_last_frame),]))^2))*0.46)
  }
}

######## Cell_Line_x_TGF #########
nobject = lengths(dataplot_Cell_Line_x_TGF)
px_distance_Cell_Line_x_TGF = NULL

for(i in 1:length(nobject)) {
  px_distance_Cell_Line_x_TGF[[i]] = data.frame(object = rep(NA,nobject[i]), distance = rep(NA,nobject[i]), 
                                       replicate = rep(NA,nobject[i]), speed = rep(NA,nobject[i]), 
                                       area = rep(NA,nobject[i]), perimeter = rep(NA,nobject[i]),
                                       displacement = rep(NA,nobject[i]))
  for (j in 1:nobject[i]) {
    if(nrow(dataplot_Cell_Line_x_TGF[[i]][[j]]) < 2) next
    ix_last_frame = nrow(dataplot_Cell_Line_x_TGF[[i]][[j]])
    px_distance_Cell_Line_x_TGF[[i]][j,] = c(j, (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_TGF[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46, i, 
                                    (sum(sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_TGF[[i]][[j]] [,c("m.cx","m.cy")]))^2))))*0.46 / (nrow(dataplot_Cell_Line_x_TGF[[i]][[j]]) / 6),
                                    mean(dataplot_Cell_Line_x_TGF[[i]][[j]][,"s.area"]), mean(dataplot_Cell_Line_x_TGF[[i]][[j]][,"s.perimeter"]),
                                    sqrt(rowSums(diff(as.matrix(dataplot_Cell_Line_x_TGF[[i]][[j]] [,c("m.cx","m.cy")][c(1,ix_last_frame),]))^2))*0.46)
  }
}

################################### Save Results ##########################################
save(px_distance_Cell_Line_x_ctrl, px_distance_Cell_Line_x_EGF, px_distance_Cell_Line_x_TGF,
     file = paste0("./RData/", outputfilename))
############################################################################################


if(file.exists(paste0(workdir,"/RData/",outputfilename))) {
cat(paste0("✅ Script completed successfully!\nResults saved to './RData/",outputfilename,"'.\n"))
} else cat("❌ Script failed. Please check if input files are provided correctly.\n")
