##################################################################
## This R script reads in the seismic traces and envelopes and  ##
## computes the covariates from the seismic traces.             ## ##################################################################

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

freq_cov_names <- paste("freq",cov_names, sep = "_")

cep_cov_names <- paste("cep",cov_names, sep = "_")

full_cov_names <- c(cov_names, freq_cov_names, cep_cov_names)

cov_df <- matrix(NA, nrow = length(decim_id), ncol = length(full_cov_names))

temp.time <- proc.time()[3]

for(j in 1:length(decim_id)){
  
  i <- decim_id[j]
  
  # Temporal domain:
  
  subset_data <- meta_df[((i-window_length+1):i), ]
  
  temp_readings <- subset_data$Data
  
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
  
  
  # Frequency domain:   
  
  temp_ts <- ts(temp_readings, frequency = 100)
  
  temp_spectrum <- spectrum(temp_ts, na.action = na.remove, plot = FALSE) 
  # Will not be so accurate if have too many NAs; tried lomb-scargle periodogram but too time-consuming.
  # For na.remove.ts this changes the “intrinsic” time scale. It is assumed that both, the new and the old time scale are synchronized at the first and the last valid observation. In between, the new series is equally spaced in the new time scale.
  temp_freq <- temp_spectrum$freq
  temp_spec <- temp_spectrum$spec
  freq_length <- length(temp_spec)
  
  temp_freq_mean <- mean(temp_spec)
  temp_freq_sd <- sd(temp_spec) 
  temp_freq_skew <- skewness(temp_spec)
  temp_freq_kurtosis <- kurtosis(temp_spec)
  temp_freq_energy <- sum(temp_spec^2)
  temp_freq_ioce <- sum((temp_spec^2)*(temp_freq)/temp_freq_energy)
  temp_freq_rms_bw <- sqrt(sum((temp_spec^2)*((temp_freq)^2)/temp_freq_energy) - temp_freq_ioce^2)
  suppressWarnings(
    temp_freq_mean_skew <- sqrt(sum((((temp_freq)-temp_freq_ioce)^3)*(temp_spec^2)/(temp_freq_energy*(temp_freq_rms_bw^3)))) # Can be NaN because taking sqrt of negative number.
  )
  temp_freq_mean_kurtosis <- sqrt(sum((((temp_freq)-temp_freq_ioce)^4)*(temp_spec^2)/(temp_freq_energy*(temp_freq_rms_bw^4)))) 
  suppressWarnings(
    temp_freq_shannon <- entropy(temp_spec, unit = "log2")
  )
  temp_freq_roa <- max(diff(temp_spec))/window_length
  temp_freq_max <- max(temp_spec)
  temp_freq_min <- min(temp_spec)
  
  # Cepstral domain:
  
  temp_ceps <- ceps(na.remove(temp_ts), qlim = c(0, 10), plot = FALSE) # Will not be so accurate if have too many NAs; similar problem to fft/spectrum computation with NAs.
  temp_quef <- temp_ceps[, 1]
  temp_amp <- temp_ceps[, 2]
  cep_length <- length(temp_amp)
  
  temp_cep_mean <- mean(temp_amp)
  temp_cep_sd <- sd(temp_amp) 
  temp_cep_skew <- skewness(temp_amp)
  temp_cep_kurtosis <- kurtosis(temp_amp)
  temp_cep_energy <- sum(temp_amp^2)
  temp_cep_ioce <- sum((temp_amp^2)*(temp_quef)/temp_cep_energy)
  temp_cep_rms_bw <- sqrt(sum((temp_amp^2)*((temp_quef)^2)/temp_cep_energy) - temp_cep_ioce^2)
  suppressWarnings(
    temp_cep_mean_skew <- sqrt(sum((((temp_quef)-temp_cep_ioce)^3)*(temp_amp^2)/(temp_cep_energy*(temp_cep_rms_bw^3))))
  )# Can be NaN because taking sqrt of negative number.
  temp_cep_mean_kurtosis <- sqrt(sum((((temp_quef)-temp_cep_ioce)^4)*(temp_amp^2)/(temp_cep_energy*(temp_cep_rms_bw^4)))) # Can be NaN because taking sqrt of negative number.
  suppressWarnings(
    temp_cep_shannon <- entropy(temp_amp, unit = "log2")
  )
  temp_cep_roa <- max(diff(temp_amp))/window_length
  temp_cep_max <- max(temp_amp)
  temp_cep_min <- min(temp_amp)
  
  cov_df[j, ] <- c(temp_mean, temp_sd, temp_skew, temp_kurtosis, temp_energy, temp_ioce, temp_rms_bw, temp_mean_skew, temp_mean_kurtosis, temp_shannon, temp_roa, temp_max, temp_min, temp_freq_mean, temp_freq_sd, temp_freq_skew, temp_freq_kurtosis, temp_freq_energy, temp_freq_ioce, temp_freq_rms_bw, temp_freq_mean_skew, temp_freq_mean_kurtosis, temp_freq_shannon, temp_freq_roa, temp_freq_max, temp_freq_min, temp_cep_mean, temp_cep_sd, temp_cep_skew, temp_cep_kurtosis, temp_cep_energy, temp_cep_ioce, temp_cep_rms_bw, temp_cep_mean_skew, temp_cep_mean_kurtosis, temp_cep_shannon, temp_cep_roa, temp_cep_max, temp_cep_min)
  
  if(j%%1000 == 0){
    print(paste(j, "/", length(decim_id), " completed.", sep = ""))
  }
  
}

time.taken <- proc.time()[3] - temp.time #2s for each 10s reading.
# 2h 20 min for lag 1h.

cov_df <- as.data.frame(cov_df)
colnames(cov_df) <- full_cov_names

# Add ratio-based covariates.

cov_df$ratio_mom <- cov_df$max/cov_df$mean
cov_df$freq_ratio_mom <- cov_df$freq_max/cov_df$freq_mean
cov_df$cep_ratio_mom <- cov_df$cep_max/cov_df$cep_mean

# Save model_df for 0.1-20 Hz only:

# unlagged_model_df <- cbind(model_df, cov_df)
# save(unlagged_model_df, file = "D:/Documents/Imperial_NTU_collaboration/Seismic data/Piton/unlagged_model_df_0120.RData")

# Save cov_df for 0.1-1, 1-5, 5-15, 0.01-.20Hz:

save(cov_df, file = "Y:/non_events/NEvent_3/unlagged_model_df_hp001.RData")

