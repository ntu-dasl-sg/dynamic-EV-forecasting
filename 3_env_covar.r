##################################################################
## This R script computes the covariates from the envelopes and ##
## combines these with those computed from the seismic traces   ##
## to form the full datasets for fitting the dynamic extreme    ##
## value model.                                                 ##
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

## 1. Load the eruption index (envelope) and the seismic data. 

# =========== Choose freq band here:
load(file = "Y:/non_events/NEvent_3/eruption_index_hp001.RData") 

miniseed <- TRUE # TRUE = other than 01-2010 event.

if(miniseed){
  
  # =========== Choose freq band here:
  load(file = "Y:/non_events/NEvent_3/NEvent_Piton_052010_traces.RData")
  
}else{
  
  # Read in unfiltered data for header:
  uv05_raw <- read.table(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/UV05_unfiltered.txt", header = TRUE)
  
  # Read in frequency-filtered data to save as trace class:  
  # ======= Change filepath accordingly (repeat for different frequency bands):
  uv05_unfilt <- read.table(file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/UV05_1-0.1Hz_new.txt", header = TRUE) 
  
  colnames(uv05_unfilt) <- colnames(uv05_raw)
  
  uv05_unfilt$DateTime <- as.POSIXct(uv05_unfilt$Time, origin="1970-01-01", tz = "GMT")
  
}

## 2. Set decimination resolution (here 1 min resolution).

if(miniseed){
  
  # =========== Choose freq band here:
  meta_df <- data.frame("Data" = tr_uv05_hp001@data)
  
}else{
  
  meta_df <- uv05_unfilt[-1, ] # Removed first observation for envelope.
  
}

meta_df$piton_env_dB <- piton_env_dB@data

decim_factor <- 100*10 # Forecast at 10 sec intervals.

# Decimination ids:

decim_id <- seq(window_length, nrow(meta_df), by = decim_factor)

model_df <- meta_df[decim_id, ]

## 3. Construct covariates.

# We construct covariates using moving windows of 1 hour aligned to the right.

window_length <- 100*60*60 # Based on past 1 hour 100Hz data to compute covariates.
# If we use the past 30 min 100Hz data to compute covariates, use 100*60*30 instead.

cov_names <- c("mean", "sd", "skew", "kurtosis", "energy", "ioce", "rms_bw", "mean_skew", "mean_kurtosis", "shannon", "roa", "max", "min")

cov_df <- matrix(NA, nrow = length(decim_id), ncol = length(cov_names))

temp.time <- proc.time()[3]

for(j in 1:length(decim_id)){
  
  i <- decim_id[j]
  
  # Temporal domain:
  
  subset_data <- meta_df[((i-window_length+1):i), ]
  temp_readings <- subset_data$piton_env_dB # Compute covariates on envelope.
  temp_mean <- mean(temp_readings, na.rm = TRUE)
  temp_sd <- sd(temp_readings, na.rm = TRUE)
  temp_skew <- skewness(temp_readings, na.rm = TRUE)
  temp_kurtosis <- kurtosis(temp_readings, na.rm = TRUE)
  temp_energy <- sum(temp_readings^2, na.rm = TRUE)
  temp_ioce <- sum((temp_readings^2)*(1:window_length)/temp_energy, na.rm = TRUE)
  temp_rms_bw <- sqrt(sum((temp_readings^2)*((1:window_length)^2)/temp_energy, na.rm = TRUE) - temp_ioce^2)
  suppressWarnings(
    temp_mean_skew <- sqrt(sum((((1:window_length)-temp_ioce)^3)*(temp_readings^2)/(temp_energy*(temp_rms_bw^3)))) 
  )
  # Can be NaN because taking sqrt of negative number.
  temp_mean_kurtosis <- sqrt(sum((((1:window_length)-temp_ioce)^4)*(temp_readings^2)/(temp_energy*(temp_rms_bw^4)), na.rm = TRUE)) 
  suppressWarnings(
    temp_shannon <- entropy(na.omit(temp_readings), unit = "log2")
  )
  temp_roa <- max(diff(temp_readings), na.rm = TRUE)/window_length
  temp_max <- max(temp_readings, na.rm = TRUE)
  temp_min <- min(temp_readings, na.rm = TRUE)
  
  # temp_cep_min <- min(temp_amp)
  
  cov_df[j, ] <- c(temp_mean, temp_sd, temp_skew, temp_kurtosis, temp_energy, temp_ioce, temp_rms_bw, temp_mean_skew, temp_mean_kurtosis, temp_shannon, temp_roa, temp_max, temp_min)

  if(j%%1000 == 0){
    print(paste(j, "/", length(decim_id), " completed.", sep = ""))
  }
  
}

time.taken <- proc.time()[3] - temp.time #About 58 min for window 1 h, 38 min for w30.

cov_df <- as.data.frame(cov_df)
colnames(cov_df) <- full_cov_names

# Add ratio-based covariates.

cov_df$env_ratio_mom <- cov_df$env_max/cov_df$env_mean # Compute on envelope.

# Save cov_df:

save(cov_df, file = "Y:/non_events/NEvent_3/unlagged_model_df_011_env.RData")


## 4. Combine covariates computed from the trace and envelope (same covariate window: 1 hour).

# Set forecast horizon by lagging covariates:
# Lag (in terms of 10 sec interval).
lag <- 60*6 # 1 hour.

# Read in all 1 hour window covariates:

load(paste("Y:/non_events/NEvent_2/unlagged_model_df_011.RData", sep = ""))

cov_df_011 <- cov_df 

colnames(cov_df_011) <- paste(colnames(cov_df_011), "_011", sep = "")

load(paste("Y:/non_events/NEvent_2/unlagged_model_df_15.RData", sep = ""))

cov_df_15 <- cov_df 

colnames(cov_df_15) <- paste(colnames(cov_df_15), "_15", sep = "")

load(paste("Y:/non_events/NEvent_2/unlagged_model_df_0120.RData", sep = ""))

cov_df_0120 <- cov_df 

colnames(cov_df_0120) <- paste(colnames(cov_df_0120), "_0120", sep = "")

load(paste("Y:/non_events/NEvent_2/unlagged_model_df_515.RData", sep = ""))

cov_df_515 <- cov_df 

colnames(cov_df_515) <- paste(colnames(cov_df_515), "_515", sep = "")

load(paste("Y:/non_events/NEvent_2/unlagged_model_df_hp001.RData", sep = ""))

cov_df_hp001 <- cov_df 

colnames(cov_df_hp001) <- paste(colnames(cov_df_hp001), "_hp001", sep = "")

# Read in all 1 hour envelope covariates

load("Y:/non_events/NEvent_2/unlagged_model_df_011_env.RData")

cov_df_011_env <- cov_df 

colnames(cov_df_011_env) <- paste(colnames(cov_df_011_env), "_011", sep = "")

load("Y:/non_events/NEvent_2/unlagged_model_df_15_env.RData")

cov_df_15_env <- cov_df 

colnames(cov_df_15_env) <- paste(colnames(cov_df_15_env), "_15", sep = "")

load("Y:/non_events/NEvent_2/unlagged_model_df_0120_env.RData")

cov_df_0120_env <- cov_df 

colnames(cov_df_0120_env) <- paste(colnames(cov_df_0120_env), "_0120", sep = "")

load("Y:/non_events/NEvent_2/unlagged_model_df_515_env.RData")

cov_df_515_env <- cov_df 

colnames(cov_df_515_env) <- paste(colnames(cov_df_515_env), "_515", sep = "")

load("Y:/non_events/NEvent_2/unlagged_model_df_hp001_env.RData")

cov_df_hp001_env <- cov_df 

colnames(cov_df_hp001_env) <- paste(colnames(cov_df_hp001_env), "_hp001", sep = "")

cov_df <- cbind(cov_df_011, cov_df_15, cov_df_0120, cov_df_515, cov_df_hp001,
                cov_df_011_env, cov_df_15_env, cov_df_0120_env, cov_df_515_env,
                cov_df_hp001_env)

lag1h_model_df <- cbind(model_df[(1+lag):nrow(model_df), ], cov_df[1:(nrow(model_df)-lag), ]) 

head(lag1h_model_df)

save(lag1h_model_df, file = "Y:/non_events/NEvent_2/lag1h_model_df_011.RData")