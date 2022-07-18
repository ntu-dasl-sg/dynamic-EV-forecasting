##################################################################
## This R script forecasts for test events/non-events using the ##
## fitted dynamic extreme value model.                          ##
################################################################## 

library(evd)
library(ggplot2)
library(MASS)
library(gridExtra)
library(car)
library(xtable)
library(generalhoslem) # Hosmer-Lemeshow Test
library(IRISSeismic)

## 1. Read in test data (to evaluate model).

# ================ Change frequency band here:
load(file = "Y:/test_data/lag1h_model_df_15.RData")
# For non-events
# load(file = "Y:/non_events/NEvent_3/lag1h_model_df_515.RData")

test_df <- lag1h_model_df

## 2. Set the index threshold.

u_chosen <- 85

#Times of exceedance in test data:
exceed_id <- which(test_df$piton_env_dB >= u_chosen)
length(exceed_id)
test_df$exceed <- 0
test_df$exceed[exceed_id] <- 1

## 3. Implement the covariate transformations.

cov_col <- which(!(colnames(test_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB", "exceed")))
length(cov_col) # 280.
cov_names <- colnames(test_df)[cov_col]

# Load previous transformation choices:

# ================ Change frequency band here:
load(file = "Y:/training_data/trans_df_15.RData") 
head(trans_df)

# Check order of names:
cov_names == trans_df$covariate

temp_df <- test_df

for (i in 1:length(cov_names)){ 
  
  covariate <- temp_df[, cov_names[i]]
  
  if(trans_df$bc_trans[i] == 1){ # Do Box-Cox transformation.
    
    temp_cov <- covariate
    
    lambda <- trans_df$lambda[i]
    buffer <- trans_df$buffer[i]
    
    if(lambda!= 1 & !is.na(lambda)){
      
      if(lambda == 0 & !is.na(lambda) & buffer == 0){
        temp_cov <- log(covariate)
      }
      if(lambda != 0 & !is.na(lambda) & buffer == 0){
        temp_cov <- (covariate)^lambda
      }
      covariate <- temp_cov
      
    }
    
  }
  
  # Do standardisation:
  
  temp_cov <- (covariate - trans_df$mean[i])/trans_df$sd[i]      
  
  temp_df[, cov_names[i]] <- temp_cov
  
}

test_df <- temp_df

# colMeans(test_df[, cov_names])

apply(test_df[, cov_names], MARGIN = 2, sd)

## 4. Forecast threshold exceedances with fitted logistic regression.

te_test_df <- test_df[, which(!(colnames(test_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB")))]

# Read in previously fitted model:
# ============== Change frequency band here:

load("Y:/training_data/te_model_15_lbound.RData")

summary(te_model)

te_predict <- predict(te_model, newdata = te_test_df, type = "response")

hist(te_predict)

te_test_df$te_prob <- te_predict

## 5. Plot forecast probabilities (select test_start and test_end based on test event/non-event).

# For 10/14/2010 event:
test_start <- as.POSIXct("2010-10-13 02:00:00", tz = "GMT")
test_end <- as.POSIXct("2010-10-15 00:00:00", tz = "GMT")
eruption_onset_start <- as.POSIXct("2010-10-14 15:20:00", tz = "GMT")
seismic_crisis_start <- as.POSIXct("2010-10-14 09:45:00", tz = "GMT")
seismic_swarm_start <- as.POSIXct("2010-10-14 10:50:00", tz = "GMT")
seismic_swarm_end <- as.POSIXct("2010-10-14 11:30:00", tz = "GMT")

# For 11/30/2009 non-event:
# test_start <- as.POSIXct("2009-11-30 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2009-12-01 00:00:00", tz = "GMT")

# For 12/22/2009 non-event:
# test_start <- as.POSIXct("2009-12-22 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2009-12-24 00:00:00", tz = "GMT")

# For 05/08/2010 non-event:
# test_start <- as.POSIXct("2009-05-08 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2009-05-10 00:00:00", tz = "GMT")

test_period <- seq(test_start, test_end, by = 10) # 10 sec decimination.
te_test_df$DateTime <- test_period

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/te_prob_15_test.png", width = 1800, height = 1200, res = 300)

plot(te_test_df$DateTime - 60*60, te_test_df$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date", xaxt = "n", ylim = c(0, 1))
# axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3] - 60*60,
# labels = c("May 08", "May 09", "May 10"))

# For 10/14/2010 event:
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3] - 60*60,
labels = c("Oct 13", "Oct 14", "Oct 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/db_envelope_test.png", width = 1800, height = 1200, res = 300)

plot(te_test_df$DateTime, lag1h_model_df$piton_env_dB, type = 'l', ylab = "Envelope (Decibels)", xlab = "Date", xaxt = "n", ylim = c(15, 100))
abline(h = 85, lty = 2, col = 2)
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)
axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3],
     labels = c("Oct 13", "Oct 14", "Oct 15"))

dev.off()

## 6. Explore major covariates (based on coefficients of the logistic regression).

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_test.png", width = 1800, height = 1200, res = 300)

plot(te_test_df$DateTime - 60*60, te_test_df$cep_kurtosis_0120, type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Date", xaxt = "n", ylim = c(-8, 4))

# For 10/14/2010 event:
axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3] - 60*60,
     labels = c("Oct 13", "Oct 14", "Oct 15"))

dev.off()

# Zoomed in:

start_buffer <- 6*60*6 # 6 hours before the eruption onset.
end_buffer <- 6*60*6 # 6 hours after the eruption onset.

plot_period <- which(te_test_df$DateTime==eruption_onset_start)

plot_period <- c((plot_period-start_buffer+1):(plot_period+end_buffer-1))

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_kurtosis_0120_test_zoom.png", width = 1800, height = 1200, res = 300)

plot(te_test_df[plot_period, "DateTime"] - 60*60, te_test_df$cep_kurtosis_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral kurtosis", xlab = "Time (UTC) on Oct 14", ylim = c(-8, 4))

# For 10/14/2010 event:
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_test.png", width = 1800, height = 1200, res = 300)

plot(te_test_df$DateTime - 60*60, te_test_df$cep_skew_0120, type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Date", xaxt = "n", ylim = c(-8, 4))

# For 10/14/2010 event:
axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3] - 60*60,
     labels = c("Oct 13", "Oct 14", "Oct 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/cep_skewness_0120_test_zoom.png", width = 1800, height = 1200, res = 300)

plot(te_test_df[plot_period, "DateTime"] - 60*60, te_test_df$cep_skew_0120[plot_period], type = 'l', ylab = "0.1-20Hz cepstral skewness", xlab = "Time (UTC) on Oct 14", ylim = c(-8, 4))

# For 10/14/2010 event:
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)

dev.off()


png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_test.png", width = 1800, height = 1200, res = 300)

plot(te_test_df$DateTime - 60*60, te_test_df$energy_hp001, type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Date", xaxt = "n", ylim = c(-2, 3))

axis(1, at = te_test_df[, "DateTime"][seq(1, nrow(te_test_df), length = 3)][1:3] - 60*60,
     labels = c("Oct 13", "Oct 14", "Oct 15"))

dev.off()

png(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/energy_hp001_test_zoom.png", width = 1800, height = 1200, res = 300)

plot(te_test_df[plot_period, "DateTime"] - 60*60, te_test_df$energy_hp001[plot_period], type = 'l', ylab = "High pass 0.01Hz energy", xlab = "Time (UTC) on Oct 14", ylim = c(-2, 3))

# For 10/14/2010 event:
abline(v = eruption_onset_start, col = 'blue')
abline(v = seismic_crisis_start, col = 'blue', lty = 3)
abline(v = seismic_swarm_start, col = 'blue', lty = 2)
abline(v = seismic_swarm_end, col = 'blue', lty = 2)

dev.off()
