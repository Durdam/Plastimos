
# Description: Use Plasticity_Index equation parameters in tSNE and visualize cells in lower dimension
# to see how the cells are clustered and weather or not the parameters are able to distinguish between
# epithelial and mesenchymal cells across various conditions.

########################################################################
########### Clustering based on Plasticity score/parameters ############
########################################################################

library("xlsx")
library("ggplot2")
library("Rtsne")

cat("Current working directory is set to:", getwd(), "\n")
# Set working directory to the root of the project
workdir = getwd()
setwd(workdir)

load("./RData/EMTresults.RData")

featNames = c("CL_GrowthRate","meanEccentricity", "meanCompactness", "Displacement_Distance",
              "AvgSpeed", "corr_Ecc_SpeedPf_Spearman")

metaDataNames = c("Plasticity_v01", "Cell_Line", "RepObj", "CL_RepObj", "Displacement", 
                  "TotalDistance", "Top20q_byAvgSpeed", "Plasticity_State_v01")

metadata = EMTcell[, metaDataNames]
dataFeat = EMTcell[, featNames]
ClusterData = scale(as.matrix(dataFeat))

cat("Setting Seed for reproducibility of results...\n")

set.seed(206)  # For reproducibility

################### tsne ################
#########################################
cat("Starting t-SNE algorithnm on Plasticity_Index parameters....\n")

# Calculate tSNE using Rtsne(0 function) 
tsne_out_01 = Rtsne(ClusterData, perplexity = 50, normalize = FALSE, verbose = TRUE, num_threads = 0)

cat("t-SNE algorithnm finished.....\n")

# Conversion of matrix to dataframe
tsne_plot_01 = data.frame(x = tsne_out_01$Y[,1], y = tsne_out_01$Y[,2])

# Plot the tsne result with metadata info
tsne_plot_01 = cbind(tsne_plot_01, metadata)
tsne_plot_01$Cell_Line = as.factor(tsne_plot_01$Cell_Line)
tsne_plot_01$Mesenchymal_State = "Low Mesenchymal"
tsne_plot_01$Mesenchymal_State[grep("High Mesenchymal", tsne_plot_01$Plasticity_State_v01)] = "High Mesenchymal"

color_CL = c("darkseagreen1", "limegreen", "darkgreen")
color_MesenchymalState = c("darkgreen", "darkgoldenrod")
color_Plasticity_State_v01 = c("steelblue3", "skyblue", "grey80")

cat("Plotting results in './images' directory.....\n")

p = ggplot(tsne_plot_01, aes(x = x, y = y, color = Cell_Line)) +
  geom_point() + scale_color_manual("Cell Lines", values = color_CL) +
  theme_classic() + 
  labs(title = paste0("tSNE: All Cells")) + xlab("tSNE_1") + ylab("tSNE_2") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
# p
ggsave(filename = paste0(workdir,"/images/tSNE_2D_AllCell_color_by_CellLine.png"),
       p, height = 6, width = 10, units = "in", dpi = 300)


p = ggplot(tsne_plot_01, aes(x = x, y = y, color = Mesenchymal_State)) +
  geom_point() + scale_color_manual("Mesenchymal State", values = color_MesenchymalState) +
  theme_classic() + 
  labs(title = paste0("tSNE: Plasticity Mesenchymal State - All Cells")) + xlab("tSNE_1") + ylab("tSNE_2") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
# p
ggsave(filename = paste0(workdir,"/images/tSNE_2D_AllCell_perplexity_50_MesenchymalState.png"),
       p, height = 6, width = 10, units = "in", dpi = 300)


p = ggplot(tsne_plot_01, aes(x = x, y = y, color = Plasticity_State_v01)) +
  geom_point() + scale_color_manual("Plasticity State", values = color_Plasticity_State_v01) +
  theme_classic() + 
  labs(title = paste0("tSNE: Plasticity State - All Cells")) + xlab("tSNE_1") + ylab("tSNE_2") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
# p
ggsave(filename = paste0(workdir,"/images/tSNE_2D_AllCell_perplexity_50_Plasticity_State.png"),
       p, height = 6, width = 12, units = "in", dpi = 300)

###############################################################################

cat("Running PCA.....\n")
dfpca = ClusterData
colnames(dfpca) = c("Growth_Rate", "meanEccentricity", "meanCompactness", "Displacement/Distance", "Average_Speed", "corr_Eccentricity_SpeedPerFrame")
res.pca = prcomp(dfpca, scale = FALSE)
fviz_eig(res.pca)
png(file = paste0(workdir,"/images/PCA_variableContribution.png"), width = 7, height = 6, res = 300, units = "in")

p = fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
print(p)
dev.off()

cat("Saving results.....\n")

save(tsne_plot_01, tsne_out_01, metadata, dataFeat, ClusterData, res.pca, 
     file = "./RData/ClusterData.RData")

write.xlsx(tsne_plot_01, file = "./output/tSNE_Clustering_Data.xlsx",
           sheetName = "tSNE_Clustering_Result", append = FALSE, row.names = FALSE)


############### script finished successfully ###################
cat(paste0("✅ Script completed successfully!\ntSNE output saved to './RData/ClusterData.RData' \noutput saved to './output/tSNE_Clustering_Data.xlsx' and \nplots saved to './images'.\n"))
################################################################