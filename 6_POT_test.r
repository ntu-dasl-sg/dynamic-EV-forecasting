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
load(file = "Y:/non_events/NEvent_3/lag1h_model_df_515.RData")

test_df <- lag1h_model_df

## 2. Set the index threshold.

u_chosen <- 82

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
load(file = "Y:/training_data/trans_df_515.RData") 
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

load("Y:/training_data/te_model_515_lbound.RData")

summary(te_model)

te_predict <- predict(te_model, newdata = te_test_df, type = "response")

hist(te_predict)

te_test_df$te_prob <- te_predict

## 5. Plot forecast probabilities (select test_start and test_end based on test event/non-event).

# For 10/14/2010 event:
# test_start <- as.POSIXct("2010-10-13 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2010-10-15 00:00:00", tz = "GMT")
# eruption_onset_start <- as.POSIXct("2010-10-14 15:20:00", tz = "GMT")
# seismic_crisis_start <- as.POSIXct("2010-10-14 09:45:00", tz = "GMT")
# seismic_swarm_start <- as.POSIXct("2010-10-14 10:50:00", tz = "GMT")
# seismic_swarm_end <- as.POSIXct("2009-10-14 11:30:00", tz = "GMT")

# For 11/30/2009 non-event:
# test_start <- as.POSIXct("2009-11-30 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2009-12-01 00:00:00", tz = "GMT")

# For 12/22/2009 non-event:
# test_start <- as.POSIXct("2009-12-22 02:00:00", tz = "GMT")
# test_end <- as.POSIXct("2009-12-24 00:00:00", tz = "GMT")

# For 05/08/2010 non-event:
test_start <- as.POSIXct("2009-05-08 02:00:00", tz = "GMT")
test_end <- as.POSIXct("2009-05-10 00:00:00", tz = "GMT")

test_period <- seq(test_start, test_end, by = 10) # 10 sec decimination.
te_test_df$DateTime <- test_period

plot(te_test_df$DateTime - 60*60, te_test_df$te_prob, type = 'l', ylab = "1 hour ahead forecast probability", xlab = "Date")
# For 10/14/2010 event:
# abline(v = eruption_onset_start, col = 'blue')
# abline(v = seismic_crisis_start, col = 'blue', lty = 3)
# abline(v = seismic_swarm_start, col = 'blue', lty = 2)
# abline(v = seismic_swarm_end, col = 'blue', lty = 2)
