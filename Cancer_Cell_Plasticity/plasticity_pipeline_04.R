
# Description: Here we measure the Plasticity_Index for each each cell using all the features we calculated so far
# in the previous pipeline scripts.
# We calculated two version of the Plasticity_Index:
# Version_01: All 6 parameters have equal weight (emtval_v01)
# Version_02: Some features like (Displacement/Distance) and Average_speed were assigned additional weight coefficients (emtval_v02)
# For the final version of this work we use the version_01.
# Further we take a 20% cutoff (+/-) to classify cells in various plasticity state using the Plasticity_Index.
# This cutoff was taken by observing the imaging data.
# Visualize the final readouts of the pipeline and compare the cell lines and various treated conditions

########################################################################
############## Compute Plasticity/EMT Index per cell ###################
########################################################################

library("ggplot2")
library("factoextra")
library("xlsx")

cat("Current working directory is set to:", getwd(), "\n")
# Set working directory to the root of the project
workdir = getwd()
setwd(workdir)

load("./RData/GrowthRate.RData")
load("./RData/All_Processed_Cell_Line_x.RData")

Plasticity_State_1 = "Low Mesenchymal / High Epithelial"
Plasticity_State_2 = "High Mesenchymal - Morphology Dependent"
Plasticity_State_3 = "High Mesenchymal - Morphology Independent"

list_AllCell = split(dfAll, f = paste0(dfAll$RepObj,"_",dfAll$Cell_Line))
EMTcell = data.frame(cellname = names(list_AllCell), Plasticity_v01 = NA, Plasticity_v02 = NA, CL_GrowthRate = NA,
                     meanEccentricity = NA, Perimeter_Area = NA, Displacement_Distance = NA, weighted_Displacement_Distance = NA,
                     AvgSpeed = NA, weighted_AvgSpeed = NA, corr_Ecc_SpeedPf_Spearman = NA)

Mes_Morph_Dep_threshold = c(80,85,90,95,100)
Mes_Morph_InDep_threshold = c(0,5,10,15,20)

####################################################################

cat("Computing Plasticity_Index for each cell.....\n.....\n")

####################################################################
for (i in 1:length(list_AllCell)) {
  
  df = list_AllCell[[i]]
  cellname = names(list_AllCell)[i]
  GR = GrowthRate_model_Output$SplineMeanGrowthFit[which(GrowthRate_model_Output$Cell_Line == unique(df$Cell_Line))]
  Eccentricity = mean(df$m.eccentricity)
  Perimeter_Area = mean(df$s.perimeter / sqrt(df$s.area))
  Displacement_Distance = unique(df$Displacement / df$TotalDistance)
  AvgSpeed = unique(df$AvgSpeed)
  SpCor_Ecc_Speed = unique(df$corr_Ecc_SpeedPf_spearman)
  
  ### Weights for the parameters
  w_Ecc = 1
  w_PA = 1
  w_DD = 2
  w_AS = 3
  w_CorES = 1
  w_GR = 1
  
  # version 01 of the Plasticity equation
  emtval_v01 = Eccentricity * Perimeter_Area * Displacement_Distance * AvgSpeed * SpCor_Ecc_Speed * GR
  
  # weighted version ofthe Plasticity equation
  weighted_Displacement_Distance = w_DD * Displacement_Distance
  weighted_AvgSpeed = w_AS * AvgSpeed
  
  emtval_v02 = (w_Ecc * Eccentricity) * (w_PA * Perimeter_Area) * 
    (w_DD * Displacement_Distance) * (w_AS * AvgSpeed) * 
    (w_CorES * SpCor_Ecc_Speed) * (w_GR * GR)
  
  # Save the output of the equation and all the varibles used to calculate the Plasticity value
  EMTcell$cellname[i] = cellname
  EMTcell$Plasticity_v01[i] = emtval_v01
  EMTcell$Plasticity_v02[i] = emtval_v02
  EMTcell$CL_GrowthRate[i] = GR
  EMTcell$meanEccentricity[i] = Eccentricity
  EMTcell$Perimeter_Area[i] = Perimeter_Area
  EMTcell$Displacement_Distance[i] = Displacement_Distance
  EMTcell$weighted_Displacement_Distance[i] = weighted_Displacement_Distance
  EMTcell$AvgSpeed[i] = AvgSpeed
  EMTcell$weighted_AvgSpeed[i] = weighted_AvgSpeed
  EMTcell$corr_Ecc_SpeedPf_Spearman[i] = SpCor_Ecc_Speed
  
}

# head(EMTcell)

CLmake = data.frame(do.call('rbind', strsplit(as.character(EMTcell$cellname),'_',fixed=TRUE)))
ncolrm = ncol(CLmake)
Cell_Line = apply(CLmake[,3:ncolrm], 1, function(row) paste(row, collapse = "_"))
RepObj = apply(CLmake[,1:2], 1, function(row) paste(row, collapse = "_"))
EMTcell$Cell_Line = Cell_Line
EMTcell$RepObj = RepObj

EMTcell$CL_RepObj = paste0(EMTcell$Cell_Line,"_",EMTcell$RepObj)
ix1 = which(EMTcell$CL_RepObj %in% dfAll$CL_RepObj)
ix2 = match(EMTcell$CL_RepObj, dfAll$CL_RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(EMTcell$CL_RepObj[ix1] == dfAll$CL_RepObj[ix3]))
# print("Perfect match in same order")

EMTcell$Displacement[ix1] = dfAll$Displacement[ix3]
EMTcell$TotalDistance[ix1] = dfAll$TotalDistance[ix3]
EMTcell$Top20q_byAvgSpeed[ix1] = dfAll$Top20q_byAvgSpeed[ix3]

### Plasticity Cutoff ####
PlasticityCutoff = data.frame(quantile(EMTcell$Plasticity_v01, probs = seq(0,1,0.05)), 
                              quantile(EMTcell$Plasticity_v02, probs = seq(0,1,0.05)))
colnames(PlasticityCutoff) = c("Plasticity_Score_v01", "Plasticity_Score_v02")
PlasticityCutoff$Quantile = rownames(PlasticityCutoff)
rownames(PlasticityCutoff) = NULL
ix_positive_mesenchymal = which(as.numeric(gsub("%", "", PlasticityCutoff$Quantile)) %in% Mes_Morph_Dep_threshold)
ix_negative_mesenchymal = which(as.numeric(gsub("%", "", PlasticityCutoff$Quantile)) %in% Mes_Morph_InDep_threshold)
## mesenchymal cells cutoff: top and bottom 10 % cells
PlasticityCutoff$Plasticity_State_v01 = Plasticity_State_1
PlasticityCutoff$Plasticity_State_v01[ix_positive_mesenchymal] = Plasticity_State_2
PlasticityCutoff$Plasticity_State_v01[ix_negative_mesenchymal] = Plasticity_State_3

Threshold_positive = min(PlasticityCutoff$Plasticity_Score_v01[which(PlasticityCutoff$Plasticity_State_v01 == Plasticity_State_2)])
Threshold_negative = max(PlasticityCutoff$Plasticity_Score_v01[which(PlasticityCutoff$Plasticity_State_v01 == Plasticity_State_3)])

### Add plasticity State to EMTcell data ####
EMTcell$Plasticity_State_v01 = Plasticity_State_1
EMTcell$Plasticity_State_v01[which(EMTcell$Plasticity_v01 >= Threshold_positive)] = Plasticity_State_2
EMTcell$Plasticity_State_v01[which(EMTcell$Plasticity_v01 <= Threshold_negative)] = Plasticity_State_3

################################################################
cat("Plotting Plasticity_Index results....\n")

# Plasticity version 01
p = ggplot(EMTcell, aes(x = Cell_Line, y = Plasticity_v01, fill = Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = c(Threshold_positive, Threshold_negative), linetype = "dashed", color = "red") + 
  labs(title = paste0("Plasticity - All Cells")) + xlab("Cell Line") + ylab("Plasticity") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
# p

ggsave(filename = paste0(workdir,"/images/Boxplot_Plasticity_AllCells.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# Plasticity version 02
# p = ggplot(EMTcell, aes(x = Cell_Line, y = Plasticity_v02, fill = Cell_Line)) +
#   geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
#   theme(legend.position = "none") +
#   labs(title = paste0("Plasticity weighted- All Cells")) + xlab("Cell Line") + ylab("Plasticity") +
#   theme(axis.text=element_text(size=12),
#         axis.title=element_text(size=16),
#         plot.title = element_text(size=18))
# p

# ggsave(filename = paste0(workdir,"images/Boxplot_weighted_Plasticity_AllCells.png"),
#        p, height = 6, width = 9, units = "in", dpi = 300)

# Generate a density plot for EMT scores
p = ggplot(EMTcell, aes(x = Plasticity_v01, color = Cell_Line)) +
  geom_density(size = 1) + xlim(-3,3) +
  theme_bw() + scale_fill_brewer(palette="Dark2") +
  labs(title = "Density Plot of EMT Scores", x = "EMT Score", y = "Density") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
# p
ggsave(filename = paste0(workdir,"/images/EMTscore_Density_Plot.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

### Plasticity State v01 All Cells ###
p = ggplot(EMTcell, aes(x=Cell_Line, y=Plasticity_v01, fill=Plasticity_State_v01)) + 
  geom_boxplot(position=position_dodge(1)) + ylim(-10,10)
p = p + theme_bw() + scale_fill_manual(values=c("lightblue3", "steelblue3","orange")) +
  labs(title = paste0("Plasticity of all cells states")) + xlab("Compound") + ylab("Plasticity") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.position="bottom")
# p
ggsave(filename = paste0(workdir,"/images/Boxplot_PlasticityStates_AllCells.png"),
       p, height = 7, width = 12, units = "in", dpi = 300)


Fraction_Plasticity = data.frame(table(paste0(EMTcell$Cell_Line,"|",EMTcell$Plasticity_State_v01)))
Fraction_Plasticity = cbind(data.frame(do.call('rbind', strsplit(as.character(Fraction_Plasticity$Var1),'|',fixed=TRUE))), Fraction_Plasticity[,2])
colnames(Fraction_Plasticity) = c("Cell_Line", "Plasticity_State_v01", "Number_of_Cells")
TotalCells_Plasticity = data.frame(table(EMTcell$Cell_Line))
colnames(TotalCells_Plasticity) = c("Cell_Line", "Number_of_Cells")
ix1 = which(Fraction_Plasticity$Cell_Line %in% TotalCells_Plasticity$Cell_Line)
ix2 = match(Fraction_Plasticity$Cell_Line, TotalCells_Plasticity$Cell_Line)
ix3 = ix2[!is.na(ix2)]
all(Fraction_Plasticity$Cell_Line[ix1] == TotalCells_Plasticity$Cell_Line[ix3])
Fraction_Plasticity$TotalCells[ix1] = TotalCells_Plasticity$Number_of_Cells[ix3]
Fraction_Plasticity$FractionCells = (Fraction_Plasticity$Number_of_Cells/Fraction_Plasticity$TotalCells) * 100
Fraction_Plasticity$Mesenchymal_State = "Low Mesenchymal"
Fraction_Plasticity$Mesenchymal_State[grep("High Mesenchymal", Fraction_Plasticity$Plasticity_State_v01)] = "High Mesenchymal"

uniCL = unique(Fraction_Plasticity$Cell_Line)
nuniCL = length(uniCL)
list_Fraction_Plasticity = split(Fraction_Plasticity, f = Fraction_Plasticity$Cell_Line)
PlasticityState = data.frame(Cell_Line = uniCL, Plasticity_State_Score = NA)
for (i in 1:nuniCL) {
  
  df = list_Fraction_Plasticity[[i]]
  cl = uniCL[i]
  
  Total_HighMesenchymal = sum(df$Number_of_Cells[which(df$Mesenchymal_State == "High Mesenchymal")])
  Total_Cells = unique(df$TotalCells)
  PlasticityState$Plasticity_State_Score[i] = Total_HighMesenchymal/Total_Cells
  PlasticityState$Plasticity_State_Score_Percent[i] = (Total_HighMesenchymal/Total_Cells) * 100
}


#### Plotting of the fraction of High-Low Mesenchymal cells in all cell lines #####
p = ggplot(PlasticityState, aes(x=Cell_Line, y=Plasticity_State_Score, size=Plasticity_State_Score)) +
  geom_point() + ylim(0,1) +
  theme_bw() +
  labs(title = paste0("Cell Line Plasticity Score")) + xlab("Cell Line") + ylab("Plasticity Score") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.position="none")
# p

ggsave(filename = paste0(workdir,"/images/Cell_Line_Plasticity_Score.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

### plot for the fraction of plastic cells ###
p = ggplot(data=Fraction_Plasticity, aes(x=Cell_Line, y=signif(FractionCells,3), fill=Plasticity_State_v01)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_text(aes(label=signif(FractionCells,3)), vjust=-0.3,
            position = position_dodge(0.9), size=3.5) +
  scale_fill_brewer(palette="Paired")+
  theme_bw() +
  labs(title = paste0("Cell Line Plasticity Fraction")) + xlab("Cell Line") + ylab("Plasticity Fraction of Cells") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.position="bottom",
        legend.text=element_text(size=12)) +
  labs(colour = NULL)
# p
ggsave(filename = paste0(workdir,"/images/Cell_Line_Plasticity_Fraction_v01.png"),
       p, height = 6, width = 14, units = "in", dpi = 300)

#### Add cell state Information to the all Data #######
ix1 = which(dfAll$CL_RepObj %in% EMTcell$CL_RepObj)
ix2 = match(dfAll$CL_RepObj, EMTcell$CL_RepObj)
ix3 = ix2[!is.na(ix2)]
if(all(dfAll$CL_RepObj[ix1] == EMTcell$CL_RepObj[ix3])) 
# print("Perfect match in same order")
dfAll$Plasticity_State_v01[ix1] = EMTcell$Plasticity_State_v01[ix3]

#############################################
############## Save results #################
#############################################

save(EMTcell, PlasticityCutoff, Threshold_negative, Threshold_positive, dfAll,
     Fraction_Plasticity, PlasticityState, 
     file = "./RData/EMTresults.RData")

write.xlsx(EMTcell, file = "./output/EMTcell.xlsx",
           sheetName = "EMTcell", append = FALSE, row.names = FALSE)

write.xlsx(PlasticityCutoff, file = "./output/PlasticityCutoff.xlsx",
           sheetName = "PlasticityCutoff", append = FALSE, row.names = FALSE)

write.table(dfAll, file = "./output/All_Data_PostProcessed.tsv", col.names = TRUE, row.names = FALSE, sep = "\t")

write.xlsx(Fraction_Plasticity, file = "./output/Plasticity_Fraction.xlsx",
           sheetName = "Plasticity_Fraction", append = FALSE, row.names = FALSE)

write.xlsx(PlasticityState, file = "./output/Plasticity_Score.xlsx",
           sheetName = "Plasticity_Score", append = FALSE, row.names = FALSE)


############### script finished successfully ###################
cat(paste0("✅ Script completed successfully!\nResults saved to './RData/EMTresults.RData' --> Input to plasticity_pipeline_05.R \nOutput saved to: './output/EMTcell.xlsx', './output/PlasticityCutoff.xlsx', './output/All_Data_PostProcessed.tsv', './output/Plasticity_Fraction.xlsx', './output/Plasticity_Score.xlsx' and Plots saved to: './images'.\n"))
################################################################