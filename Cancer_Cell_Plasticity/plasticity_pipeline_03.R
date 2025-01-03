
# Description: We tested 2 models to fit and measure growth rate for the various conditions.
# We tested a linear model and what worked best for us was a spline function.
# We used the spline function to finally get growth rate parameter (G) that we used in our Plasticity_Index equation

########################################################################
######################### Compute growth rate ##########################
########################################################################

library("xlsx")
library("ggplot2")

cat("Current working directory is set to:", getwd(), "\n")
# Set working directory to the root of the project
workdir = getwd()
setwd(workdir)

load("./RData/All_Processed_Cell_Line_x.RData")

######################## Argument inputs #############################
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Parse named arguments

args_list = data.frame(do.call('rbind', strsplit(as.character(args),'=',fixed=TRUE)))
args_list$X1 = trimws(args_list$X1)
args_list$X2 = trimws(args_list$X2)

nframe = as.integer(args_list$X2[which(args_list$X1 == "nframe")])
nframe_gap = as.integer(args_list$X2[which(args_list$X1 == "nframe_gap")])
time_0_frame = as.character(args_list$X2[which(args_list$X1 == "time_0_frame")])

# print(nframe)
# print(length(nframe))

# Check if 'nframe' and 'nframe_gap' are provided
if (length(nframe == 1) && length(nframe_gap == 1) && length(time_0_frame == 1)) {
  cat("Input data check: ✅  'nframe', 'nframe_gap' and 'time_0_frame' arguments are required and provided.\n")
} else {
  cat("❌  Error: 'nframe', 'nframe_gap' and 'time_0_frame' arguments are required. Use nframe=<value> nframe_gap=<value> time_0_frame=<value>.\n")
}

# Print the arguments
cat("📊 Number of frames (nframe):", nframe, "\n")
cat("⏳ Gap between frames (nframe_gap):", nframe_gap, "\n")
cat("⏲️ Time = 0 frame included (time_0_frame):", time_0_frame, "\n")

# Compute total time based on whether time = 0 frame is included
if (time_0_frame == "y") {
  total_time <- (nframe - 1) * nframe_gap
} else {
  total_time <- nframe * nframe_gap
}
total_time = total_time/60

# Print the result
cat("🕒 Total time:", total_time, "hours\n\n")

##########################################
# Calcuate Growth rate for each Cell line

df_growthrate = data.frame(table(paste0(dfAll$Cell_Line,"|",dfAll$frame)))
# head(df_growthrate)
df_growthrate = cbind(data.frame(do.call('rbind', strsplit(as.character(df_growthrate$Var1),'|',fixed=TRUE))), df_growthrate[,2])
colnames(df_growthrate) = c("Cell_Line", "Frame", "Number_of_Cells")
df_growthrate$Frame = as.numeric(df_growthrate$Frame)

ixframe_above_151 = which(df_growthrate$Frame > 151)
if(length(ixframe_above_151) > 0) { df_growthrate = df_growthrate[-ixframe_above_151, ] }

df_growthrate = df_growthrate[order(df_growthrate$Cell_Line, df_growthrate$Frame), ]
nCellLine = length(unique(df_growthrate$Cell_Line))
df_growthrate$Time =  rep(seq(0, total_time*60, by = nframe_gap) / 60, nCellLine)

# Plot the cell counts over time
p = ggplot(df_growthrate, aes(x = Time, y = Number_of_Cells, color = Cell_Line)) +
  geom_line() +
  labs(title = "Cell Growth Over Time", x = "Time (hours)", y = "Number of Cells") +
  scale_fill_brewer(palette="Dark2") + theme_bw() + xlim(0,total_time) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
# p
ggsave(filename = paste0(workdir,"/images/Cell_Growth_over_time.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

cl1 = "Cell_Line_x_ctrl"
cl2 = "Cell_Line_x_EGF"
cl3 = "Cell_Line_x_TGF"

ixcl1 = which(df_growthrate$Cell_Line==cl1)
ixcl2 = which(df_growthrate$Cell_Line==cl2)
ixcl3 = which(df_growthrate$Cell_Line==cl3)

# Log-transform the cell counts for all cell lines
log_cl1 = log(df_growthrate$Number_of_Cells[ixcl1])
log_cl2 = log(df_growthrate$Number_of_Cells[ixcl2])
log_cl3 = log(df_growthrate$Number_of_Cells[ixcl3])

df_growthrate$log_nCells = c(log_cl1, log_cl2, log_cl3)
Time = df_growthrate$Time[ixcl1]

#######################################################
################# linear model ########################
#######################################################

print("Computing Growth Rate: fitting linear model.....")

# Fit linear models (on the log-transformed data)
fit_cl1 = lm(log_cl1 ~ Time)
fit_cl2 = lm(log_cl2 ~ Time)
fit_cl3 = lm(log_cl3 ~ Time)

# Extract the growth rates (slope of the fitted line)
growth_rate_cl1_exp = coef(fit_cl1)[2]
growth_rate_cl2_exp = coef(fit_cl2)[2]
growth_rate_cl3_exp = coef(fit_cl3)[2]

## save the single value growth rates in a data.frame
GrowthRate_model_Output = data.frame(Cell_Line = c(cl1, cl2, cl3),
                                     ExponentialGrowthfit = c(growth_rate_cl1_exp, growth_rate_cl2_exp, growth_rate_cl3_exp))

# You can visualize the fitted lines on log-transformed data
p = ggplot(df_growthrate, aes(x = Time)) +
  geom_point(aes(y = log_nCells, color = Cell_Line), alpha = 0.5) +
  geom_smooth(aes(y = log_nCells, color = Cell_Line), method = "lm", se = FALSE) +
  labs(title = "Exponential Growth Model", x = "Time (hours)", y = "Log(Number of Cells)") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))


ggsave(filename = paste0(workdir,"/images/Exponential_Cell_Growth_Fitted.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

cat("Linear model plot saved to: './images/Exponential_Cell_Growth_Fitted.png'\n")

cells_line1 = df_growthrate$Number_of_Cells[ixcl1]
cells_line2 = df_growthrate$Number_of_Cells[ixcl2]
cells_line3 = df_growthrate$Number_of_Cells[ixcl3]

###############################################################
######################## Split fitting ########################
###############################################################

print("Computing Growth Rate: fitting Spline function.....")

# Fit a smooth spline to the cell data for both lines
spline_fit_line1 <- smooth.spline(Time, cells_line1, spar = 0.8)
spline_fit_line2 <- smooth.spline(Time, cells_line2, spar = 0.8)
spline_fit_line3 <- smooth.spline(Time, cells_line3, spar = 0.8)

# Predict fitted values
fitted_spline_line1 <- predict(spline_fit_line1, Time)$y
fitted_spline_line2 <- predict(spline_fit_line2, Time)$y
fitted_spline_line3 <- predict(spline_fit_line3, Time)$y

df_growthrate$Spline_fit = c(fitted_spline_line1, fitted_spline_line2, fitted_spline_line3)

# Differentiate the spline to get the rate of growth (slope) at each time point
spline_deriv_line1 <- predict(spline_fit_line1, Time, deriv = 1)$y
spline_deriv_line2 <- predict(spline_fit_line2, Time, deriv = 1)$y
spline_deriv_line3 <- predict(spline_fit_line3, Time, deriv = 1)$y

# Calculate doubling time for logistic growth
growth_rate_line1 = mean(spline_deriv_line1 / (log(2) * fitted_spline_line1))
growth_rate_line2 = mean(spline_deriv_line2 / (log(2) * fitted_spline_line2))
growth_rate_line3 = mean(spline_deriv_line3 / (log(2) * fitted_spline_line3))

GrowthRate_model_Output$SplineMeanGrowthFit = c(growth_rate_line1, growth_rate_line2, growth_rate_line3)

# visualize the fitted lines on spline growth rate data
color_CL = c("darkseagreen1", "limegreen", "darkgreen")
p = ggplot(df_growthrate, aes(x = Time)) +
  geom_point(aes(y = Number_of_Cells, color = Cell_Line), alpha = 0.75) +
  geom_line(aes(y = Spline_fit, color = Cell_Line), size = 1) +
  labs(x = "Time (hours)", y = "Number of Cells") +
  theme_bw() +scale_color_manual("Cell Lines", values = color_CL) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=10)))

# p
ggsave(filename = paste0(workdir,"/images/SplineFit_model_Cell_Growth.png"),
       p, height = 6, width = 10, units = "in", dpi = 300)

cat("Spline model plot saved to: './images/SplineFit_model_Cell_Growth.pdf'.\n")

################ Save Results #####################

save(df_growthrate, GrowthRate_model_Output,
     cl1, cl2, cl3,
     file = "./RData/GrowthRate.RData")

write.xlsx(GrowthRate_model_Output, file = "./output/GrowthRate_model_Output.xlsx",
           sheetName = "GrowthRate_perCellLine", append = FALSE, row.names = FALSE)

write.xlsx(df_growthrate, file = "./output/GrowthRate_fitdata.xlsx",
           sheetName = "GrowthRate_fitdata", append = FALSE, row.names = FALSE)

############### script finished successfully ###################
cat(paste0("✅ Script completed successfully!\nResults saved to './RData/GrowthRate.RData' \noutput saved to './output/GrowthRate_model_Output.xlsx' and './output/GrowthRate_fitdata.xlsx'\nplots saved to './images'.\n"))

