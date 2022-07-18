##################################################################
## This R script fits the dynamic extreme value model to        ##
## training data using step-wise variable selection.            ##
################################################################## 

library(evd)
library(ggplot2)
library(MASS)
library(gridExtra)
library(car)
library(xtable)
library(generalhoslem) # Hosmer-Lemeshow Test

## 1. Read in training data.

# ================ Choose frequency band here:

load(file = "Y:/training_data/Train_Piton_112009/lag1h_model_df_15.RData") 

lag1h_model_df_train1 <- lag1h_model_df

load(file = "Y:/training_data/Train_Piton_122009/lag1h_model_df_15.RData") 

lag1h_model_df_train2 <- lag1h_model_df

load(file = "Y:/training_data/Train_Piton_012010/lag1h_model_df_15.RData") 

lag1h_model_df_train3 <- lag1h_model_df

chosen_col <- colnames(lag1h_model_df_train1)

# ================ Choose number of training sets here:
model_df <- rbind(lag1h_model_df_train1, lag1h_model_df_train2[, chosen_col], lag1h_model_df_train3[, chosen_col])

# model_df <- lag1h_model_df_train3[, chosen_col]

train1_end <- nrow(lag1h_model_df_train1)
train2_end <- train1_end + nrow(lag1h_model_df_train2)
train3_end <- nrow(model_df)


## 2. Set an eruption index and threshold.

eruption_index <- model_df$piton_env_dB
u_chosen <- 85

# Times of exceedance
exceed_id <- which(eruption_index >= u_chosen)
length(exceed_id)
model_df$exceed <- 0
model_df$exceed[exceed_id] <- 1

## 3. Transform covariates via Box-Cox transforms and standardise them.

cov_col <- which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB", "exceed")))

length(cov_col) # 280.
cov_names <- colnames(model_df)[cov_col]

# Box-Cox transformation for more normal distribution shapes:

lambda_vector <- rep(NA, length(cov_col))
buffer_vector <- rep(NA, length(cov_col))

test_col <- c(1:length(cov_col)) # Skip any columns if required.


for (i in test_col){ 
  
  covariate <- model_df[, cov_col[i]]
  
  if(sum(!is.nan(covariate))>0){
    
    buffer <- -min(covariate, na.rm = TRUE)
    if(buffer>0){
      bc_result <- boxcox(covariate + 1.1*buffer ~ 1, plot = FALSE)
      buffer_vector[i] <- buffer
    }else{
      
      if(buffer == 0){
        buffer <- sort(covariate[covariate!=0], decreasing = FALSE)[1]
        bc_result <- boxcox(covariate + 1.1*buffer ~ 1, plot = FALSE)
        buffer_vector[i] <- buffer
      }else{bc_result <- boxcox(covariate ~ 1, plot = FALSE)}
    }
    # Round lambda to the nearest 0.5.
    bc_lambda <- round(bc_result$x[bc_result$y == max(bc_result$y)]*2)/2
    
    lambda_vector[i] <- bc_lambda    
    
  }
  
}

# Only transform positive variables to avoid buffer choice:

buffer_vector[is.na(buffer_vector)] <- 0

for (i in 1:length(cov_col)){ 
  
  covariate <- model_df[, cov_col[i]]
  
  lambda <- lambda_vector[i]
  buffer <- buffer_vector[i]
  
  if(lambda!= 1 & !is.na(lambda)){
    
    temp_cov <- covariate
    
    if(lambda == 0 & !is.na(lambda) & buffer == 0){
      temp_cov <- log(covariate)
    }
    if(lambda != 0 & !is.na(lambda) & buffer == 0){
      temp_cov <- (covariate)^lambda
    }
    model_df[, cov_col[i]] <- temp_cov
    
  }
  
}

# Save transformation choices:
trans_df <- data.frame("covariate" = cov_names, "lambda" = lambda_vector, "buffer" = buffer_vector)
trans_df$bc_trans <- 0
trans_df$bc_trans[!is.na(trans_df$lambda) & trans_df$lambda!=1 & trans_df$buffer == 0] <- 1

# Standardize to common scale:

mean_vector <- rep(NA, length(cov_col))
sd_vector <- rep(NA, length(cov_col))

for (i in 1:length(cov_col)){ 
  
  covariate <- model_df[, cov_col[i]]
  
  temp_mean <- mean(covariate, na.rm = TRUE)
  temp_sd <- sd(covariate, na.rm = TRUE)
  mean_vector[i] <- temp_mean
  sd_vector[i] <- temp_sd
  
  stand_cov <- (covariate - temp_mean)/temp_sd
  
  model_df[, cov_col[i]] <- stand_cov
  
}

trans_df$mean <- mean_vector
trans_df$sd <- sd_vector

head(trans_df)

# ================ Choose frequency band here:
save(trans_df, file = "Y:/training_data/trans_df_15.RData")

# Check for NaNs and remove covariates with NaNs: 
nan_col <- which(lapply(model_df, FUN = function(x){sum(is.nan(x))}) > 0)
model_df <- model_df[, !(colnames(model_df) %in% names(nan_col))]

## 4. Fit logistic regression for threshold exceedances

te_df <- model_df[, which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB")))]

# Check for multicollinearity in covariates:
# Order covariates in terms of how much they inform univariate models.

exceed_vec <- te_df$exceed
aic_vec <- rep(NA, (ncol(te_df)-1))

for(i in 1:(ncol(te_df)-1)){
  temp_cov <- te_df[, i]
  te_univariate <- glm(exceed_vec ~ temp_cov, family = binomial(link = "logit"))
  aic_vec[i] <- te_univariate$aic  
}

cov_order <- order(aic_vec, decreasing = FALSE)
cor_mat <- cor(as.matrix(te_df[, 1:(ncol(te_df)-1)]))

# Remove covariates with 0.6 or higher absolute correlation with others:

cov_remove <- 0
cov_keep <- 0
cov_omit <- 0

for(i in 1:length(cov_order)){
  
  if(!(cov_order[i] %in% cov_remove)){
    
    cov_row <- cor_mat[cov_order[i], ]
    if(i>1){
      cov_omit <- cov_order[1]:cov_order[i-1] 
    }
    cov_exceed <- which(abs(cov_row) > 0.6)
    cov_remove <- c(cov_remove, cov_exceed[cov_exceed!=cov_order[i] & !(cov_exceed %in% cov_omit)])    
    
  }    
}

cov_remove <- unique(cov_remove)
cov_remove <- cov_remove[-1]

cov_remove

cor_mat[-cov_remove, -cov_remove]

te_df <- te_df[, -cov_remove]

head(te_df)

te_logistic <- glm(exceed ~ ., family = binomial, data = te_df)

summary(te_logistic)

# Conduct step-wise variable selection:

te_stepAIC <- stepAIC(te_logistic, direction = "both")

te_model <- te_stepAIC

te_summary <- summary(te_model)

# ================ Choose frequency band here:

# Obtain LaTeX summary table.
print(xtable(te_summary, type = "latex"), file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/training_data/te_summary_15_lbound.tex")
# Save model summary.
save(te_model, file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/training_data/te_model_15_lbound.RData")

# If saved earlier, 
load(file = "Y:/training_data/te_model_15_lbound.RData")
# For training event 3 only: 
# load(file = "Y:/archive/index_trace_env_15/te_model.RData")

## 5. Compute and plot training forecasts against event timings.

te_predict <- predict(te_model, type = "response")

hist(te_predict)

model_df$te_prob <- te_predict

# ================ Choose number of training sets here:
forecast_1 <- model_df[1:train1_end, ]
forecast_2 <- model_df[(train1_end+1):train2_end, ]
forecast_3 <- model_df[(train2_end+1):train3_end, ]

##################
# Training set 1 #
##################

# For 11/2009 event:
eruption_onset_start <- as.POSIXct("2009-11-05 17:00:00", tz = "GMT")
seismic_crisis_start <- as.POSIXct("2009-11-05 15:30:00", tz = "GMT")
seismic_swarm_start <- as.POSIXct("2009-11-05 15:40:00", tz = "GMT")
seismic_swarm_end <- as.POSIXct("2009-11-05 16:10:00", tz = "GMT")

# One hour for covariate window; one hour lag.
time_start <- as.POSIXct("2009-11-04 02:00:00", tz = "GMT")
time_end <- as.POSIXct("2009-11-07 00:00:00", tz = "GMT")

time_period <- seq(time_start, time_end, by = 10) # 10 sec decimination.
length(time_period) #25561

nrow(lag1h_model_df_train1) # 25201

lag1h_model_df_train1$DateTime <- time_period

# Plot of full event period:

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train1_lbound.png", width = 1800, height = 1200, res = 300)

# Note: Forecast for 02:00 made at 01:00.
plot(lag1h_model_df_train1[, "DateTime"] - 60*60, forecast_1$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date", xaxt = "n", ylim = c(0, 1))
axis(1, at = lag1h_model_df_train1[, "DateTime"][seq(1, nrow(lag1h_model_df_train1), length = 4)][1:3] - 60*60, 
     labels = c("Nov 04", "Nov 05", "Nov 06"))

dev.off()

# Plot zoomed into day of eruption:

start_buffer <- 6*60*3 # 3 hours before the eruption onset.
end_buffer <- 6*60*3 # 3 hours after the eruption onset.

plot_period <- which(lag1h_model_df_train1$DateTime==eruption_onset_start)

plot_period <- c((plot_period-start_buffer+1):(plot_period+end_buffer-1))

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train1_zoom_lbound.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train1[plot_period, "DateTime"] - 60*60, forecast_1$te_prob[plot_period], type = 'l',  ylab = "1 hour ahead forecast probability", xlab = "Time (UTC) on Nov 05", ylim = c(0, 1))
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
dev.off()

# Major covariates:

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train1.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train1$DateTime - 60*60, forecast_1$cep_kurtosis_0120, type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Date", xaxt = "n", ylim = c(-8, 4)) 
axis(1, at = lag1h_model_df_train1[, "DateTime"][seq(1, nrow(lag1h_model_df_train1), length = 4)][1:3] - 60*60, 
     labels = c("Nov 04", "Nov 05", "Nov 06"))

dev.off()

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train1_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train1[plot_period, "DateTime"] - 60*60, forecast_1$cep_kurtosis_0120[plot_period], type = 'l',  ylab = "0.1-20Hz cepstral kurtosis", xlab = "Time (UTC) on Nov 05", ylim = c(-8, 4)) 
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
dev.off()


png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train1.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train1$DateTime - 60*60, forecast_1$cep_skew_0120, type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Date", xaxt = "n", ylim = c(-8, 4))
axis(1, at = lag1h_model_df_train1[, "DateTime"][seq(1, nrow(lag1h_model_df_train1), length = 4)][1:3] - 60*60, 
     labels = c("Nov 04", "Nov 05", "Nov 06"))

dev.off()

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train1_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train1[plot_period, "DateTime"] - 60*60, forecast_1$cep_skew_0120[plot_period], type = 'l',  ylab = "0.1-20Hz cepstral skewness", xlab = "Time (UTC) on Nov 05", ylim = c(-8, 4)) 
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train1.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train1$DateTime - 60*60, forecast_1$energy_hp001, type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Date", xaxt = "n", ylim = c(-2, 3))
axis(1, at = lag1h_model_df_train1[, "DateTime"][seq(1, nrow(lag1h_model_df_train1), length = 4)][1:3] - 60*60, 
     labels = c("Nov 04", "Nov 05", "Nov 06"))

dev.off()

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train1_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train1[plot_period, "DateTime"] - 60*60, forecast_1$energy_hp001[plot_period], type = 'l',  ylab = "High pass 0.01Hz energy", xlab = "Time (UTC) on Nov 05", ylim = c(-2, 3))
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
dev.off()

##################
# Training set 2 #
##################

# For 12/2009 event:
eruption_onset_start <- as.POSIXct("2009-12-14 14:40:00", tz = "GMT")
seismic_crisis_start <- as.POSIXct("2009-12-14 13:30:00", tz = "GMT")
seismic_swarm_start <- as.POSIXct("2009-12-14 13:40:00", tz = "GMT")
seismic_swarm_end <- as.POSIXct("2009-12-14 14:12:00", tz = "GMT")

time_start <- as.POSIXct("2009-12-13 02:00:00", tz = "GMT")
time_end <- as.POSIXct("2009-12-15 00:00:00", tz = "GMT")

time_period <- seq(time_start, time_end, by = 10) # 10 sec decimination.
length(time_period) #25561

nrow(lag1h_model_df_train2) # 25201

lag1h_model_df_train2$DateTime <- time_period

# Plot of full event period:

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train2_lbound.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2[, "DateTime"] - 60*60, forecast_2$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date", xaxt = "n", ylim = c(0, 1))
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

# Plot zoomed into day of eruption:

start_buffer <- 6*60*3 # 3 hours before the eruption onset.
end_buffer <- 6*60*3 # 3 hours after the eruption onset.

plot_period <- which(lag1h_model_df_train2$DateTime==eruption_onset_start)

plot_period <- c((plot_period-start_buffer+1):(plot_period+end_buffer-1))

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train2_zoom_lbound.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train2[plot_period, "DateTime"] - 60*60, forecast_2$te_prob[plot_period], type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Time (UTC) on Dec 14", ylim = c(0, 1))
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
dev.off()

# Major covariates:

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train2.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2$DateTime - 60*60, forecast_2$cep_kurtosis_0120, type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Date", xaxt = "n", ylim = c(-8, 4)) 
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train2_zoom.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2[plot_period, "DateTime"] - 60*60, forecast_2$cep_kurtosis_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Time (UTC) on Dec 14", ylim = c(-8, 4))

abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train2.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2$DateTime - 60*60, forecast_2$cep_skew_0120, type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Date", xaxt = "n", ylim = c(-8, 4))
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train2_zoom.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2[plot_period, "DateTime"] - 60*60, forecast_2$cep_skew_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Time (UTC) on Dec 14", ylim = c(-8, 4))

abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train2.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2$DateTime - 60*60, forecast_2$energy_hp001, type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Date", xaxt = "n", ylim = c(-2, 3)) 
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train2_zoom.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train2[plot_period, "DateTime"] - 60*60, forecast_2$energy_hp001[plot_period], type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Time (UTC) on Dec 14", ylim = c(-2, 3)) 

abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
axis(1, at = lag1h_model_df_train2[, "DateTime"][seq(1, nrow(lag1h_model_df_train2), length = 3)] - 60*60, 
     labels = c("Dec 13", "Dec 14", "Dec 15"))

dev.off()

##################
# Training set 3 #
##################

# For Jan 2010 event:
s_crisis_start <- as.POSIXct("2010-01-02 07:50:00", tz = "GMT")
s_swarm_start <- as.POSIXct("2010-01-02 08:20:00", tz = "GMT")
s_swarm_end <- as.POSIXct("2010-01-02 09:02:00", tz = "GMT")
eruption_onset_start <- as.POSIXct("2010-01-02 10:20:00", tz = "GMT")

# Plot of full event period:
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train3_lbound.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train3[, "DateTime"] - 60*60, forecast_3$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date", xaxt = "n", ylim = c(0, 1))
axis(1, at = lag1h_model_df_train3[, "DateTime"][seq(1, nrow(lag1h_model_df_train3), length = 4)][1:3] - 60*60, 
     labels = c("Jan 01", "Jan 02", "Jan 03"))

dev.off()

# Plot of full event period (for training set 3 only fit):
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_15_train3_only.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train3[, "DateTime"] - 60*60, model_df$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date", xaxt = "n", ylim = c(0, 1))
axis(1, at = lag1h_model_df_train3[, "DateTime"][seq(1, nrow(lag1h_model_df_train3), length = 4)][1:3] - 60*60,
     labels = c("Jan 01", "Jan 02", "Jan 03"))

dev.off()


# Plot zoomed into day of eruption:

start_buffer <- 6*5 # 5 minutes before the crisis.
end_buffer <- 6*60*3 # 3 hours after the eruption onset.

plot_period <- which(lag1h_model_df_train3$DateTime>=s_crisis_start & lag1h_model_df_train3$DateTime<=eruption_onset_start)

plot_period <- c((plot_period[1]-start_buffer+1):(plot_period[1]-1), plot_period, (plot_period[length(plot_period)]+ 1):(plot_period[length(plot_period)]+end_buffer-1))

# DateTime in model_df is for forecasted time so is 1 hour after time of covariates. 
png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_515_train3_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train3[plot_period, "DateTime"] - 60*60, forecast_3$te_prob[plot_period], type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Time (UTC) on Jan 02", ylim = c(0, 1))
abline(v = eruption_onset_start, col = 'blue')
abline(v = s_crisis_start, col = 'blue', lty = 3)
abline(v = s_swarm_start, col = 'blue', lty = 2)
abline(v = s_swarm_end, col = 'blue', lty = 2)
dev.off()

# For training set 3 only fit:
# png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_15_train3_only_zoom.png", width = 1800, height = 1200, res = 300)
# plot(lag1h_model_df_train3[plot_period, "DateTime"] - 60*60, model_df$te_prob[plot_period], type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Time (UTC) on Jan 02", ylim = c(0, 1))
# abline(v = eruption_onset_start, col = 'blue')
# abline(v = s_crisis_start, col = 'blue', lty = 3)
# abline(v = s_swarm_start, col = 'blue', lty = 2)
# abline(v = s_swarm_end, col = 'blue', lty = 2)
# dev.off()

# Major covariates:

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train3.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train3$DateTime - 60*60, forecast_3$cep_kurtosis_0120, type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Date", xaxt = "n", ylim = c(-8, 4)) 
axis(1, at = lag1h_model_df_train3[, "DateTime"][seq(1, nrow(lag1h_model_df_train3), length = 4)][1:3] - 60*60, 
     labels = c("Jan 01", "Jan 02", "Jan 03"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_train3_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train3[plot_period, "DateTime"] - 60*60, forecast_3$cep_kurtosis_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Time (UTC) on Jan 02", ylim = c(-8, 4)) 
abline(v = eruption_onset_start, col = 'blue')
abline(v = s_crisis_start, col = 'blue', lty = 3)
abline(v = s_swarm_start, col = 'blue', lty = 2)
abline(v = s_swarm_end, col = 'blue', lty = 2)
dev.off()


png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train3.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train3$DateTime - 60*60, forecast_3$cep_skew_0120, type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Date", xaxt = "n", ylim = c(-8, 4)) 
axis(1, at = lag1h_model_df_train3[, "DateTime"][seq(1, nrow(lag1h_model_df_train3), length = 4)][1:3] - 60*60, 
     labels = c("Jan 01", "Jan 02", "Jan 03"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_train3_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train3[plot_period, "DateTime"] - 60*60, forecast_3$cep_skew_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Time (UTC) on Jan 02", ylim = c(-8, 4)) 
abline(v = eruption_onset_start, col = 'blue')
abline(v = s_crisis_start, col = 'blue', lty = 3)
abline(v = s_swarm_start, col = 'blue', lty = 2)
abline(v = s_swarm_end, col = 'blue', lty = 2)
dev.off()


png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train3.png", width = 1800, height = 1200, res = 300)

plot(lag1h_model_df_train3$DateTime - 60*60, forecast_3$energy_hp001, type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Date", xaxt = "n", ylim = c(-2, 3)) 
axis(1, at = lag1h_model_df_train3[, "DateTime"][seq(1, nrow(lag1h_model_df_train3), length = 4)][1:3] - 60*60, 
     labels = c("Jan 01", "Jan 02", "Jan 03"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_train3_zoom.png", width = 1800, height = 1200, res = 300)
plot(lag1h_model_df_train3[plot_period, "DateTime"] - 60*60, forecast_3$energy_hp001[plot_period], type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Time (UTC) on Jan 02", ylim = c(-2, 3)) 
abline(v = eruption_onset_start, col = 'blue')
abline(v = s_crisis_start, col = 'blue', lty = 3)
abline(v = s_swarm_start, col = 'blue', lty = 2)
abline(v = s_swarm_end, col = 'blue', lty = 2)
dev.off()


## 6. Compute goodness-of-fit tests for the logistic regression.

# Deviance chi-squared test for goodness of fit:
anova(te_model, test = "Chisq")

1-pchisq(te_model$null.deviance - te_model$deviance, df=(te_model$df.null-te_model$df.residual))

# Hosmer-Lemeshow test on null of equality between expected and observed frequency of exceedance:

hl_test <- logitgof(obs = te_df$exceed, exp = model_df$te_prob, g = 10) 
# Note default number of quantiles of risk = 10 but need g = 4 to have >1 in expected frequencies table for chi-squared approximation to hold. 

hl_test
hl_test$expected

## 7. Check for conditional independence between threshold exceedances via ACF plots.

res.pearson <- residuals(te_model, type="pearson")
#Constant probability model:
cp_model <- glm(exceed ~ 1, family = binomial, data = te_df)
res.pearson_cp <- residuals(cp_model, type="pearson")

# Training set 1

acf_te_res <- acf(res.pearson[1:train1_end], plot = FALSE)
acf_cp_res <- acf(res.pearson_cp[1:train1_end], plot = FALSE)

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_te_15_train1_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cp_res, type = 'l', col = 'red', main = '', lty = 3)
lines(acf_te_res$lag, acf_te_res$acf)
legend(22.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Logistic regression", "Constant probability"), bty = "n")

dev.off()

# Training set 2

acf_te_res <- acf(res.pearson[(train1_end+1):train2_end], plot = FALSE)
acf_cp_res <- acf(res.pearson_cp[(train1_end+1):train2_end], plot = FALSE)
# barlett_ci <- qnorm(0.975)/ sqrt(nrow(te_df))

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_te_15_train2_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cp_res, type = 'l', col = 'red', main = '', lty = 3)
lines(acf_te_res$lag, acf_te_res$acf)
legend(22.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Logistic regression", "Constant probability"), bty = "n")

dev.off()

# Training set 3

acf_te_res <- acf(res.pearson[(train2_end+1):train3_end], plot = FALSE)
acf_cp_res <- acf(res.pearson_cp[(train2_end+1):train3_end], plot = FALSE)
# barlett_ci <- qnorm(0.975)/ sqrt(nrow(te_df))

# For training set 3 only fit:
# acf_te_res <- acf(res.pearson, plot = FALSE)
# acf_cp_res <- acf(res.pearson_cp, plot = FALSE)

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_te_15_train3_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cp_res, type = 'l', col = 'red', main = '', lty = 3)
lines(acf_te_res$lag, acf_te_res$acf)
legend(22.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Logistic regression", "Constant probability"), bty = "n")

dev.off()


## 8. Fit dynamic GP model to threshold excesses.

# Estimate and fix the GPD shape parameter:

test_model_1 <- fpot(model_df$piton_env_dB, threshold = u_chosen, model = "gpd", cmax = FALSE)
test_model_1
fixed_shape <- test_model_1$estimate['shape'] 

fixed_shape + c(-1, 1)*qnorm(0.975)*test_model_1$std.err['shape']

# If the estimated shape parameter is not significantly different from zero, we use the definition of the GPD distribution with shape parameter equal to 0. That is, an exponential distribution (a Gamma distribution with shape/dispersion parameter 1). Otherwise, we use a GPD regression as defined in the Bee et al. paper.

model_df$excess <- model_df$piton_env_dB - u_chosen 

excess_df <- model_df[, which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB", "te_prob", "exceed")))]

excess_df <- excess_df[excess_df$excess>=0, ]

# ================ Indicate if shape parameter is set to zero:
zero_shape <- FALSE

source('D:/Documents/Imperial_NTU_collaboration/Seismic data/dynamic-EV-forecasting/GPD_regression.r')

# Check for multicollinearity in covariates:

# Order covariates in terms of how much they inform univariate models.

exceed_vec <- excess_df$excess
aic_vec <- rep(NA, (ncol(excess_df)-1))

if(!zero_shape){
  
  start_beta <- c(5, 0)
  z_t <- excess_df$excess
  
  for(i in 1:(ncol(excess_df)-1)){
    temp_cov <-  cbind(matrix(1, nrow = nrow(excess_df), ncol = 1), as.matrix(excess_df[, i]))
    excess_univariate <- optim(start_beta, GPD_nll, X = temp_cov, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE) #, lower = lower_val, upper = upper_val)
    aic_vec[i] <- 2*excess_univariate$value + 2*length(excess_univariate$par)
  }
  
}else{
  
  for(i in 1:(ncol(excess_df)-1)){
    temp_cov <- excess_df[, i]
    excess_univariate <- glm(exceed_vec ~ temp_cov, family = Gamma(link = "log"))
    aic_vec[i] <- excess_univariate$aic  
  }
  
}

cov_order <- order(aic_vec, decreasing = FALSE)
cor_mat <- cor(as.matrix(excess_df[, 1:(ncol(excess_df)-1)]))

# Remove covariates with 0.6 or higher absolute correlation with others:

cov_remove_2 <- 0
cov_keep <- 0
cov_omit <- 0

for(i in 1:length(cov_order)){
  
  if(!(cov_order[i] %in% cov_remove_2)){
    
    cov_row <- cor_mat[cov_order[i], ]
    if(i>1){
      cov_omit <- cov_order[1]:cov_order[i-1] 
    }
    cov_exceed <- which(abs(cov_row) > 0.6)
    cov_remove_2 <- c(cov_remove_2, cov_exceed[cov_exceed!=cov_order[i] & !(cov_exceed %in% cov_omit)])    
    
  }    
}

cov_remove_2 <- unique(cov_remove_2)
cov_remove_2 <- cov_remove_2[-1]

excess_df <- excess_df[, -cov_remove_2]

# Step-wise variable selection for GPD/exponential regression.

if(zero_shape){
  
  excess_expreg <- glm(excess ~ ., family = Gamma(link = "log"), data = excess_df)
  
  summary(excess_expreg, dispersion = 1) # Dispersion parameter does not affect mean estimates but does impact their standard errors.
  
  excess_stepAIC <- stepAIC(excess_expreg, direction = "both")
  
  excess_model <- excess_stepAIC
  
  excess_AIC <- AIC(excess_model)
  
  excess_summary <- summary(excess_model, dispersion = 1)
  
  excess_pmean <- predict(excess_model, type = "response")
  
}else{
  
  temp_stepAIC <- GPD_stepAIC(excess_df, fixed_shape, max_it = 40) # GPD_stepAIC: forward then backwards; GPD_stepAIC2: backwards then forwards (as per stepAIC).
  
  excess_AIC <- min(temp_stepAIC$AIC_store)
  
  cov_id <- temp_stepAIC[['current_cov']]
  cov_id <- cov_id[1:length(temp_stepAIC[['current_cov']])] 
  cov_names_2 <- colnames(excess_df)
  temp_cov <-  cbind(matrix(1, nrow = nrow(excess_df), ncol = 1), as.matrix(excess_df[, cov_id]))
  start_beta <- c(5, rep(0, length(cov_id))) 
  excess_univariate <- optim(start_beta, GPD_nll, X = temp_cov, z_t = exceed_vec, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = TRUE) 
  
  fisher_info<-solve(excess_univariate$hessian)
  prop_sigma<-sqrt(diag(fisher_info))
  
  
  excess_summary <- data.frame(" " = c("(Intercept)", cov_names_2[cov_id]), "Estimate" = excess_univariate$par, "Std. Error" = prop_sigma)
  excess_summary$"z value" <- excess_univariate$par/prop_sigma
  excess_summary$"Pr(>|z|)" <- 1 - pnorm(abs(excess_summary$"z value"))
  
  v_t <- exp(temp_cov%*%excess_univariate$par)
  
  excess_pmean <- (v_t)/(1-fixed_shape)
  
}

# ================ Choose frequency band here:

# Obtain LaTeX summary table.
print(xtable(excess_summary, type = "latex"), file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/training_data/excess_summary_15_lbound.tex")

## 9. Check for conditional independence between threshold excesses via ACF plots:

res.excess <- excess_df$excess - excess_pmean

# Constant GP model:
cgp_model <- fpot(model_df$piton_env_dB, threshold = u_chosen, model = "gpd", cmax = FALSE, std.err = FALSE, method = "Nelder-Mead")
AIC(cgp_model) >  excess_AIC 

excess_pmean_cgp <- (cgp_model$estimate[1])/(1-cgp_model$estimate[2])
res.excess_cgp <- excess_df$excess - excess_pmean_cgp

set_id <- c(rep(1, train1_end), rep(2, train2_end - train1_end + 1), rep(3, train3_end - train2_end + 1))[model_df$excess >= 0]

# Training set 1:

acf_excess_res <- acf(res.excess[set_id == 1]) 
acf_cgp_res <- acf(res.excess_cgp[set_id == 1])

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_excess_15_train1_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cgp_res, type = 'l', col = 'red', main = '', lty = 3, xlim = c(0, 20))
lines(acf_excess_res$lag, acf_excess_res$acf)
legend(12.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Dynamic GP", "Constant GP"), bty = "n")

dev.off()

# Training set 2:

acf_excess_res <- acf(res.excess[set_id == 2]) 
acf_cgp_res <- acf(res.excess_cgp[set_id == 2])

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_excess_15_train2_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cgp_res, type = 'l', col = 'red', main = '', lty = 3, xlim = c(0, 20))
lines(acf_excess_res$lag, acf_excess_res$acf)
legend(12.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Dynamic GP", "Constant GP"), bty = "n")

dev.off()


# Training set 3:

acf_excess_res <- acf(res.excess[set_id == 3]) 
acf_cgp_res <- acf(res.excess_cgp[set_id == 3])

# For training set 3 only fit:
# acf_excess_res <- acf(res.excess) 
# acf_cgp_res <- acf(res.excess_cgp)

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/acf_excess_15_train3_lbound.png", width = 1800, height = 1200, res = 300)

plot(acf_cgp_res, type = 'l', col = 'red', main = '', lty = 3, xlim = c(0, 20))
lines(acf_excess_res$lag, acf_excess_res$acf)
legend(12.5, 0.95, col = c("black", "red"), lty = c(1, 3), legend = c("Dynamic GP", "Constant GP"), bty = "n")

dev.off()

## 10. Check goodness-of-fit of GPD regression.

# Graphical validation as in Coles (2001): Standardise excesses to standard exponentially distributed variables.

# Without assuming shape = 0, for GPD regression: 

shape_est <- fixed_shape

stand_excess <- log(1 + shape_est*(excess_df$excess/excess_pmean))/shape_est 
not_nan_id <- which(!is.nan(stand_excess))
stand_excess_ordered <- sort(stand_excess[not_nan_id], decreasing = FALSE)

# If assuming shape = 0: 
# stand_excess_2 <- excess_df$excess/excess_pmean
# stand_excess_ordered_2 <- sort(stand_excess_2, decreasing = FALSE)

exp_quantile <- -log(1 - (1:length(stand_excess))/(length(stand_excess)+1))

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/excess_qq_lbound.png", width = 1200, height = 1200, res = 300)

plot(exp_quantile[not_nan_id], stand_excess_ordered, type = 'p', ylab = "Empirical quantiles", xlab = "Theoretical quantiles", ylim = c(0, 20), xlim = c(0, 8)) #, main = "Non-zero shape parameter", ylim = c(0, 20), xlim = c(0, 20))
abline(a = 0, b = 1)

dev.off()

# Without assuming shape = 0, for constant GPD:

shape_est <- cgp_model$estimate[2]

stand_excess <- log(1 + shape_est*(excess_df$excess/excess_pmean_cgp))/shape_est 
not_nan_id <- which(!is.nan(stand_excess))
stand_excess_ordered <- sort(stand_excess[not_nan_id], decreasing = FALSE)

# If assuming shape = 0: 
# stand_excess_2 <- excess_df$excess/excess_pmean_cgp
# stand_excess_ordered_2 <- sort(stand_excess_2, decreasing = FALSE)

exp_quantile <- -log(1 - (1:length(stand_excess))/(length(stand_excess)+1))

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/excess_qq_lbound_cgpd.png", width = 1200, height = 1200, res = 300)

plot(exp_quantile[not_nan_id], stand_excess_ordered, type = 'p', ylab = "Empirical quantiles", xlab = "Theoretical quantiles", ylim = c(0, 20), xlim = c(0, 8)) #, main = "Non-zero shape parameter")
abline(a = 0, b = 1)

dev.off()

# Save forecast results.
# ================ Choose frequency band here:
save(model_df, file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/training_data/model_df_15.RData")
