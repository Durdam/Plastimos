
#####################################################
#### Plotting and visualizations of the features ####
#####################################################

library("ggplot2")
library("here")

# Set working directory to the root of the project
workdir = here()
setwd(workdir)

# make folder structures
if(!file.exists(paste0(workdir,"RData"))) dir.create(paste0(workdir,"RData"))
if(!file.exists(paste0(workdir,"output"))) dir.create(paste0(workdir,"output"))
if(!file.exists(paste0(workdir,"images"))) dir.create(paste0(workdir,"images"))

load("./RData/All_Processed_MCF7_MB231.RData")

label_AvgSpeed = paste("Average Speed ","[",intToUtf8(181),"M/hr]", sep = "")
label_TotalDistance = paste("Total Distance ","[",intToUtf8(181),"M]", sep = "")
label_SpeedPerFrame = paste("Speed Per Frame ","[",intToUtf8(181),"M/hr]", sep = "")
label_DistancePerFrame = paste("Distance Per Frame ","[",intToUtf8(181),"M]", sep = "")
label_Displacement = paste("Displacement ","[",intToUtf8(181),"M]", sep = "")

# top 20% fastest moving cells only
dftopq = dfAll[which(dfAll$Top20q_byAvgSpeed == "yes"), ]

# Eccentricity
p = ggplot(dftopq, aes(x=Cell_Line, y=m.eccentricity, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Eccentricity")) + xlab("Cell Line") + ylab("Eccentricity") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_Eccentricity.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)


# perimeter
p = ggplot(dftopq, aes(x=Cell_Line, y=s.perimeter, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Perimeter")) + xlab("Cell Line") + ylab("Perimeter") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_Perimeter.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# area
p = ggplot(dftopq, aes(x=Cell_Line, y=s.area, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Area")) + xlab("Cell Line") + ylab("Area") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_Area.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# AvgSpeed
p = ggplot(dftopq, aes(x=Cell_Line, y=AvgSpeed, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Average Speed")) + xlab("Cell Line") + ylab(label_AvgSpeed) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_AvgSpeed.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# SpeedPerFrame
p = ggplot(dftopq, aes(x=Cell_Line, y=SpeedPerFrame, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Speed Per Frame")) + xlab("Cell Line") + ylab(label_SpeedPerFrame) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_SpeedPerFrame.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# TotalDistance
p = ggplot(dftopq, aes(x=Cell_Line, y=TotalDistance, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Total Distance")) + xlab("Cell Line") + ylab(label_TotalDistance) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_TotalDistance.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# DistancePerFrame
p = ggplot(dftopq, aes(x=Cell_Line, y=DistancePerFrame, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Distance Per Frame")) + xlab("Cell Line") + ylab(label_DistancePerFrame) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_DistancePerFrame.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# Displacement
p = ggplot(dftopq, aes(x=Cell_Line, y=Displacement, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Displacement")) + xlab("Cell Line") + ylab(label_Displacement) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_Displacement.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)


# corr_Ecc_SpeedPf Spearman
p = ggplot(dftopq, aes(x=Cell_Line, y=corr_Ecc_SpeedPf_spearman, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Spearman correlation: Eccentricity-SpeedPerFrame")) + xlab("Cell Line") + ylab("r") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_corrSpearman_Ecc_SpeedPf.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

# corr_Ecc_SpeedPf Pearson
p = ggplot(dftopq, aes(x=Cell_Line, y=corr_Ecc_SpeedPf_pearson, fill=Cell_Line)) +
  geom_boxplot() + scale_fill_brewer(palette="Dark2") + theme_bw() + 
  theme(legend.position = "none") +
  labs(title = paste0("Pearson correlation: Eccentricity-SpeedPerFrame")) + xlab("Cell Line") + ylab("r") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p

ggsave(filename = paste0(workdir,"images/Boxplot_corrPearson_Ecc_SpeedPf.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)
