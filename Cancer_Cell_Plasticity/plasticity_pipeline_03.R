
# Description: We tested 3 models to fit and mesure growth rate for the various condition.
# We tested a linear model, logistic growth model and lastly what worked best for us was a spline function.
# We used the spline function to finally get growth rate parameter (G) that we used in our Plasticity_Index equation

########################################################################
######################### Compute growth rate ##########################
########################################################################

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

##########################################
# Calcuate Growth rate for each Cell line

df_growthrate = data.frame(table(paste0(dfAll$Cell_Line,"|",dfAll$frame)))
head(df_growthrate)
df_growthrate = cbind(data.frame(do.call('rbind', strsplit(as.character(df_growthrate$Var1),'|',fixed=TRUE))), df_growthrate[,2])
colnames(df_growthrate) = c("Cell_Line", "Frame", "Number_of_Cells")
df_growthrate$Frame = as.numeric(df_growthrate$Frame)
df_growthrate = df_growthrate[-which(df_growthrate$Frame > 151), ]
df_growthrate = df_growthrate[order(df_growthrate$Cell_Line, df_growthrate$Frame), ]
nCellLine = length(unique(df_growthrate$Cell_Line))
df_growthrate$Time =  rep(seq(0, 25*60, by = 10) / 60, nCellLine)

# Plot the cell counts over time
p = ggplot(df_growthrate, aes(x = Time, y = Number_of_Cells, color = Cell_Line)) +
  geom_line() +
  labs(title = "Cell Growth Over Time", x = "Time (hours)", y = "Number of Cells") +
  scale_fill_brewer(palette="Dark2") + theme_bw() + xlim(0,25) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))
p
ggsave(filename = paste0(workdir,"images/Cell_Growth_over_time.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)

cl1 = "MB231_ctrl"
cl2 = "MB231_EGF"
cl3 = "MB231_TGF"
cl4 = "MCF7_ctrl"
cl5 = "MCF7_EGF"
cl6 = "MCF7_TGF"

ixcl1 = which(df_growthrate$Cell_Line==cl1)
ixcl2 = which(df_growthrate$Cell_Line==cl2)
ixcl3 = which(df_growthrate$Cell_Line==cl3)
ixcl4 = which(df_growthrate$Cell_Line==cl4)
ixcl5 = which(df_growthrate$Cell_Line==cl5)
ixcl6 = which(df_growthrate$Cell_Line==cl6)


# Log-transform the cell counts for all cell lines
log_cl1 = log(df_growthrate$Number_of_Cells[ixcl1])
log_cl2 = log(df_growthrate$Number_of_Cells[ixcl2])
log_cl3 = log(df_growthrate$Number_of_Cells[ixcl3])
log_cl4 = log(df_growthrate$Number_of_Cells[ixcl4])
log_cl5 = log(df_growthrate$Number_of_Cells[ixcl5])
log_cl6 = log(df_growthrate$Number_of_Cells[ixcl6])

df_growthrate$log_nCells = c(log_cl1, log_cl2, log_cl3, log_cl4, log_cl5, log_cl6)
Time = df_growthrate$Time[ixcl1]

#######################################################
################# linear model ########################
#######################################################
# Fit linear models (on the log-transformed data)
fit_cl1 = lm(log_cl1 ~ Time)
fit_cl2 = lm(log_cl2 ~ Time)
fit_cl3 = lm(log_cl3 ~ Time)
fit_cl4 = lm(log_cl4 ~ Time)
fit_cl5 = lm(log_cl5 ~ Time)
fit_cl6 = lm(log_cl6 ~ Time)

# Extract the growth rates (slope of the fitted line)
growth_rate_cl1_exp = coef(fit_cl1)[2]
growth_rate_cl2_exp = coef(fit_cl2)[2]
growth_rate_cl3_exp = coef(fit_cl3)[2]
growth_rate_cl4_exp = coef(fit_cl4)[2]
growth_rate_cl5_exp = coef(fit_cl5)[2]
growth_rate_cl6_exp = coef(fit_cl6)[2]

# Print the growth rates for all cell lines
paste0("Exponential growth rate for ",cl1,": ", growth_rate_cl1_exp)
paste0("Exponential growth rate for ",cl2,": ", growth_rate_cl2_exp)
paste0("Exponential growth rate for ",cl3,": ", growth_rate_cl3_exp)
paste0("Exponential growth rate for ",cl4,": ", growth_rate_cl4_exp)
paste0("Exponential growth rate for ",cl5,": ", growth_rate_cl5_exp)
paste0("Exponential growth rate for ",cl6,": ", growth_rate_cl6_exp)

## save the single value growth rates in a data.frame
GrowthRate_model_Output = data.frame(Cell_Line = c(cl1, cl2, cl3, cl4, cl5, cl6),
                                     ExponentialGrowthfit = c(growth_rate_cl1_exp, growth_rate_cl2_exp, growth_rate_cl3_exp, growth_rate_cl4_exp, growth_rate_cl5_exp, growth_rate_cl6_exp))

# You can visualize the fitted lines on log-transformed data
p = ggplot(df_growthrate, aes(x = Time)) +
  geom_point(aes(y = log_nCells, color = Cell_Line), alpha = 0.5) +
  geom_smooth(aes(y = log_nCells, color = Cell_Line), method = "lm", se = FALSE) +
  labs(title = "Exponential Growth Model", x = "Time (hours)", y = "Log(Number of Cells)") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))


ggsave(filename = paste0(workdir,"images/Exponential_Cell_Growth_Fitted.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)


#######################################################
############# logistic growth model ###################
#######################################################
# Define a logistic growth model function
logistic_growth <- function(t, N0, K, r) {
  K / (1 + ((K - N0) / N0) * exp(-r * t))
}

# Fit the logistic model to the data for Cell Line 1
cells_line1 = df_growthrate$Number_of_Cells[ixcl1]
fit_logistic_cl1 <- nls(cells_line1 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line1), K = max(cells_line1), r = 0.1))

cells_line2 = df_growthrate$Number_of_Cells[ixcl2]
fit_logistic_cl2 <- nls(cells_line2 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line2), K = max(cells_line2), r = 0.1))

cells_line3 = df_growthrate$Number_of_Cells[ixcl3]
fit_logistic_cl3 <- nls(cells_line3 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line3), K = max(cells_line3), r = 0.1))

cells_line4 = df_growthrate$Number_of_Cells[ixcl4]
fit_logistic_cl4 <- nls(cells_line4 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line4), K = max(cells_line4), r = 0.1))

cells_line5 = df_growthrate$Number_of_Cells[ixcl5]
fit_logistic_cl5 <- nls(cells_line5 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line5), K = max(cells_line5), r = 0.1))

cells_line6 = df_growthrate$Number_of_Cells[ixcl6]
fit_logistic_cl6 <- nls(cells_line6 ~ logistic_growth(Time, N0, K, r), 
                        start = list(N0 = min(cells_line6), K = max(cells_line6), r = 0.1))

# Summary of the fit
summary(fit_logistic_cl1)
summary(fit_logistic_cl2)
summary(fit_logistic_cl3)
summary(fit_logistic_cl4)
summary(fit_logistic_cl5)
summary(fit_logistic_cl6)

# Extract fitted values
fitted_line1 <- predict(fit_logistic_cl1)
fitted_line2 <- predict(fit_logistic_cl2)
fitted_line3 <- predict(fit_logistic_cl3)
fitted_line4 <- predict(fit_logistic_cl4)
fitted_line5 <- predict(fit_logistic_cl5)
fitted_line6 <- predict(fit_logistic_cl6)

df_growthrate$LogisticGrowth_fit = c(fitted_line1, fitted_line2, fitted_line3, fitted_line4, fitted_line5, fitted_line6)

# Extract the growth rates (slope of the fitted line)
growth_rate_cl1_logistic = coef(fit_logistic_cl1)["r"]
growth_rate_cl2_logistic = coef(fit_logistic_cl2)["r"]
growth_rate_cl3_logistic = coef(fit_logistic_cl3)["r"]
growth_rate_cl4_logistic = coef(fit_logistic_cl4)["r"]
growth_rate_cl5_logistic = coef(fit_logistic_cl5)["r"]
growth_rate_cl6_logistic = coef(fit_logistic_cl6)["r"]

# Print the growth rates for all cell lines
paste0("Logistic growth rate for ",cl1,": ", growth_rate_cl1_logistic)
paste0("Logistic growth rate for ",cl2,": ", growth_rate_cl2_logistic)
paste0("Logistic growth rate for ",cl3,": ", growth_rate_cl3_logistic)
paste0("Logistic growth rate for ",cl4,": ", growth_rate_cl4_logistic)
paste0("Logistic growth rate for ",cl5,": ", growth_rate_cl5_logistic)
paste0("Logistic growth rate for ",cl6,": ", growth_rate_cl6_logistic)

GrowthRate_model_Output$LogisticGrowthFit = c(growth_rate_cl1_logistic, growth_rate_cl2_logistic, growth_rate_cl3_logistic, growth_rate_cl4_logistic, growth_rate_cl5_logistic, growth_rate_cl6_logistic)

# visualize the fitted lines on logistic growth rate data
p = p = ggplot(df_growthrate, aes(x = Time)) +
  geom_point(aes(y = Number_of_Cells, color = Cell_Line), alpha = 0.5) +
  geom_line(aes(y = LogisticGrowth_fit, color = Cell_Line), size = 1) +
  labs(title = "Cell Growth with Logistic Growth Model", x = "Time (hours)", y = "Number of Cells") +
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16),
        plot.title = element_text(size=18))

p
ggsave(filename = paste0(workdir,"images/Logistic_Growth_model_Cell_Growth_Fitted.png"),
       p, height = 6, width = 9, units = "in", dpi = 300)


###############################################################
######################## Split fitting ########################
###############################################################
# Fit a smooth spline to the cell data for both lines
spline_fit_line1 <- smooth.spline(Time, cells_line1, spar = 0.8)
spline_fit_line2 <- smooth.spline(Time, cells_line2, spar = 0.8)
spline_fit_line3 <- smooth.spline(Time, cells_line3, spar = 0.8)
spline_fit_line4 <- smooth.spline(Time, cells_line4, spar = 0.8)
spline_fit_line5 <- smooth.spline(Time, cells_line5, spar = 0.8)
spline_fit_line6 <- smooth.spline(Time, cells_line6, spar = 0.8)

# Predict fitted values
fitted_spline_line1 <- predict(spline_fit_line1, Time)$y
fitted_spline_line2 <- predict(spline_fit_line2, Time)$y
fitted_spline_line3 <- predict(spline_fit_line3, Time)$y
fitted_spline_line4 <- predict(spline_fit_line4, Time)$y
fitted_spline_line5 <- predict(spline_fit_line5, Time)$y
fitted_spline_line6 <- predict(spline_fit_line6, Time)$y

df_growthrate$Spline_fit = c(fitted_spline_line1, fitted_spline_line2, fitted_spline_line3, fitted_spline_line4, fitted_spline_line5, fitted_spline_line6)

# Differentiate the spline to get the rate of growth (slope) at each time point
spline_deriv_line1 <- predict(spline_fit_line1, Time, deriv = 1)$y
spline_deriv_line2 <- predict(spline_fit_line2, Time, deriv = 1)$y
spline_deriv_line3 <- predict(spline_fit_line3, Time, deriv = 1)$y
spline_deriv_line4 <- predict(spline_fit_line4, Time, deriv = 1)$y
spline_deriv_line5 <- predict(spline_fit_line5, Time, deriv = 1)$y
spline_deriv_line6 <- predict(spline_fit_line6, Time, deriv = 1)$y

# Calculate doubling time for logistic growth
growth_rate_line1 = mean(spline_deriv_line1 / (log(2) * fitted_spline_line1))
growth_rate_line2 = mean(spline_deriv_line2 / (log(2) * fitted_spline_line2))
growth_rate_line3 = mean(spline_deriv_line3 / (log(2) * fitted_spline_line3))
growth_rate_line4 = mean(spline_deriv_line4 / (log(2) * fitted_spline_line4))
growth_rate_line5 = mean(spline_deriv_line5 / (log(2) * fitted_spline_line5))
growth_rate_line6 = mean(spline_deriv_line6 / (log(2) * fitted_spline_line6))


# GrowthRate_model_Output$SplineMeanGrowthFit = c(mean_growth_rate_line1, mean_growth_rate_line2, mean_growth_rate_line3, mean_growth_rate_line4, mean_growth_rate_line5, mean_growth_rate_line6)
GrowthRate_model_Output$SplineMeanGrowthFit = c(growth_rate_line1, growth_rate_line2, growth_rate_line3, growth_rate_line4, growth_rate_line5, growth_rate_line6)

# visualize the fitted lines on spline growth rate data
color_CL = c("darkseagreen1", "limegreen", "darkgreen", "salmon", "plum", "darkgoldenrod")
p = ggplot(df_growthrate, aes(x = Time)) +
  geom_point(aes(y = Number_of_Cells, color = Cell_Line), alpha = 0.75) +
  geom_line(aes(y = Spline_fit, color = Cell_Line), size = 1) +
  labs(x = "Time (hours)", y = "Number of Cells") + ylim(50,450) +
  theme_bw() +scale_color_manual("Cell Lines", values = color_CL) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        plot.title = element_text(size=18),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 18)) +
  guides(colour = guide_legend(override.aes = list(size=10)))

p
ggsave(filename = paste0(workdir,"images/SplineFit_model_Cell_Growth.pdf"),
       p, height = 6, width = 10, units = "in", dpi = 300)

save(df_growthrate, GrowthRate_model_Output,
     cl1, cl2, cl3, cl4, cl5, cl6,
     file = "./RData/GrowthRate.RData")

write.xlsx(GrowthRate_model_Output, file = "./output/GrowthRate_model_Output.xlsx",
           sheetName = "GrowthRate_perCellLine", append = FALSE, row.names = FALSE)

write.xlsx(df_growthrate, file = "./output/GrowthRate_fitdata.xlsx",
           sheetName = "GrowthRate_fitdata", append = FALSE, row.names = FALSE)
