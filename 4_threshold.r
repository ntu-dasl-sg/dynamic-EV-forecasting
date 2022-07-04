##################################################################
## This R script computes the p-values and graphs relevant for  ##
## selecting a threshold of the eruption index (trace envelope) ##
## for the dynamic extreme value model.                         ##
################################################################## 

library(evd)
library(ggplot2)
library(zoo)
library(moments)
library(entropy)
library(seewave)
library(eva) # Automatic threshold selection.
library(texmex) # Extremal index for cluster detection.
library(lemon)
library(tseries)

## 1. Choose the eruption index to select a threshold for and the range of threshold values to test. 

# =========== Choose freq band here:

load(file = "Y:/training_data/Train_Piton_112009/lag1h_model_df_15.RData") 

train_env1 <- lag1h_model_df$piton_env_dB

load(file = "Y:/training_data/Train_Piton_122009/lag1h_model_df_15.RData") 

train_env2 <- lag1h_model_df$piton_env_dB

load(file = "Y:/training_data/Train_Piton_012010/lag1h_model_df_15.RData") 

train_env3 <- lag1h_model_df$piton_env_dB

# =========== Choose eruption index here:
eruption_index <- train_env3
# Use  c(train_env1, train_env2, train_env3) if computing threshold based on all training events.

# Set the range of threshold values to test:
q1 <- quantile(eruption_index, p = 0.85) 
q2 <- quantile(eruption_index, p = 0.999) 
threshold_vec <- seq(round(q1), round(q2), length.out = 25)

# Set a significance level for the goodness-of-fit tests:
sig_level <- 0.1

## 2. Conduct the Anderson-Darling (AD) test.

temp.time <- proc.time()[3]

auto_thres <- gpdSeqTests(data = eruption_index, thresholds = threshold_vec, method = "ad", nsim = 499) 

time.taken <- proc.time()[3] - temp.time 
# Training set 3: 2.2555  min; all three training sets = 18.6 hours.

auto_thres
at_df <- data.frame("Threshold" = rep(auto_thres$threshold, 2))
at_df$Variable <- NA
temp_var <- c(auto_thres$p.values, auto_thres$ForwardStop)
temp_var_name <- c(rep("p value", length(threshold_vec)), 
                   rep("ForwardStop", length(threshold_vec))
                   )
at_df$Value <- temp_var
at_df$Variable <- temp_var_name

ad_thres <- auto_thres$threshold[which(auto_thres$ForwardStop > sig_level)][1]

ad_plot <- ggplot(data = at_df) + geom_line(aes(x = Threshold, y = Value, color = Variable)) + geom_point(aes(x = Threshold, y = Value, color = Variable)) + theme_classic() + ggtitle("AD test")  + coord_cartesian(ylim = c(0, 1.5)) + geom_hline(yintercept = sig_level, lty = 2) 

## 3. Conduct the Cramer-von Mises (CVM) test.

temp.time <- proc.time()[3]

auto_thres_2 <- gpdSeqTests(data = eruption_index, thresholds = threshold_vec, method = "cvm", nsim = 499) 

# data = eruption_index_trunc or train_env1@data etc.

time.taken <- proc.time()[3] - temp.time 
# Training set 3: 2.22966 min; all three training sets = 21.53 hours.

auto_thres_2
at_df_2 <- data.frame("Threshold" = rep(auto_thres_2$threshold, 2))
at_df_2$Variable <- NA
temp_var <- c(auto_thres_2$p.values, auto_thres_2$ForwardStop)
temp_var_name <- c(rep("p value", length(threshold_vec)), 
                   rep("ForwardStop", length(threshold_vec))
)
at_df_2$Value <- temp_var
at_df_2$Variable <- temp_var_name

cvm_thres <- auto_thres_2$threshold[auto_thres_2$ForwardStop == min(auto_thres_2$ForwardStop)]

cvm_plot <- ggplot(data = at_df_2) + geom_line(aes(x = Threshold, y = Value, color = Variable)) + geom_point(aes(x = Threshold, y = Value, color = Variable)) + theme_classic() + ggtitle("CVM test")  + coord_cartesian(ylim = c(0, 1.5)) + geom_hline(yintercept = sig_level, lty = 2) 

## 4. Plot the results of both goodness-of-fit tests on the same graph.

at_df$Test <- "AD"
at_df_2$Test <- "CVM"

at_df_full <- rbind(at_df, at_df_2)

at_df_full <- at_df_full[at_df_full$Variable == "p value", ]

adcvm_plot <- ggplot(data = at_df_full) + geom_line(aes(x = Threshold, y = Value, color = Test)) + geom_point(aes(x = Threshold, y = Value, color = Test)) + theme_classic()  + coord_cartesian(ylim = c(0, 1)) + geom_hline(yintercept = sig_level, lty = 2) + ylab("p value") + theme(legend.position="bottom")

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/train3_pvalue.png", width = 1250, height = 1000, res = 300)
adcvm_plot
dev.off()

## 5. Check Generalised Pareto Distribution (GPD) goodness-of-fit plots for chosen threshold value.
at_df_full

# Fit the peak-over-threshold model with MLE optim:
# =========== Choose threshold value here:
u_chosen <- 85

test_model_2a <- fpot(eruption_index, threshold = u_chosen, model = "gpd" , cmax = FALSE, method = "Nelder-Mead") #, std.err = FALSE)
test_model_2a


par(mfrow = c(2, 2))
plot(test_model_2a)

## 6. Visualise the threshold on the eruption index and compute the number of exceedances.

# Rmd figure: fig.height = 9, fig.width = 7.

par(mfrow = c(3, 1))

plot(train_env1, type = "l", xlab = "Time index", ylab = "Training set 1: hp001", main = "dB Envelope")
abline(h = u_chosen, col='red', lwd=2)

plot(train_env2, type = "l", xlab = "Time index", ylab = "Training set 2: hp001", main = "dB Envelope")
abline(h = u_chosen, col='red', lwd=2)

plot(train_env3, type = "l", xlab = "Time index", ylab = "Training set 3: hp001", main = "dB Envelope")
abline(h = u_chosen, col='red', lwd=2)

# Number of exceedances:
sum(train_env1 > u_chosen)
sum(train_env2 > u_chosen)
sum(train_env3 > u_chosen)
