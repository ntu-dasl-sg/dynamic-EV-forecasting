##################################################################
## This R script reads in the raw seismic data and computes the ##
## trace envelope index.                                        ##
##################################################################

## 1. Load packages.

library(IRISSeismic)
library(pracma)
library(WaveletComp)

## 2. Read in data and convert to trace class to use IRISSeismic functions.

# Use demo data format for trace class:

iris <- new("IrisClient")
# Set the starttime and endtime
starttime <- as.POSIXct("2012-01-24", tz="GMT")
endtime <- as.POSIXct("2012-01-25", tz="GMT")
# Get the waveform
st <- getDataselect(iris,"AK","PIN","","BHZ",starttime,endtime)
# Get the first trace
tr1 <- st@traces[[1]]

# Read in test/training data and filter to frequency bands as required:

miniseed <- FALSE # TRUE if in miniseed format; FALSE if using pre-processed data in text file format.

if(miniseed){
  
  temp.time <- proc.time()[3]
  # ======= Change filepath accordingly:
  uv05_raw <- readMiniseedFile("Z:/non_events/NEvent_2/NEvent_Piton_122009.mseed")
  time.taken <- proc.time()[3] - temp.time # 6 minutes.
  
  # Ensure sequence has even length because FFT is very slow on odd length sequences:
  
  if(mod(length(uv05_raw@traces[[1]]), 2) == 1){
    
    # Remove first row:
    
    uv05_raw@traces[[1]]@data <- uv05_raw@traces[[1]]@data[-1]
    
  }          
  
  tr_uv05_raw <- uv05_raw@traces[[1]]  
  
  # Frequency-filter:
  
  # 1. Highpass 0.01Hz 
  
  temp_butter <- signal::butter(n = 2, W = 0.01, type = "high")
  highpass_signal <- signal::filter(temp_butter$b, temp_butter$a, tr_uv05_raw@data)
  tr_uv05_hp001 <- tr_uv05_raw
  tr_uv05_hp001@data <- as.vector(highpass_signal)
  
  # 2. 0.1-20Hz
  
  tr_uv05_0120 <- IRISSeismic::butterworth(x = tr_uv05_raw, n = 2, low = 0.1, high = 20)
  
  # 3. 1-0.1Hz
  
  tr_uv05_011 <- IRISSeismic::butterworth(x = tr_uv05_raw, n = 2, low = 0.1, high = 1)
  
  # 4. 5-15 Hz
  
  tr_uv05_515 <- IRISSeismic::butterworth(x = tr_uv05_raw, n = 2, low = 5, high = 15)   
  
  # 5. 1-5 Hz
  
  tr_uv05_15 <- IRISSeismic::butterworth(x = tr_uv05_raw, n = 2, low = 1, high = 5)   
  
  # ======= Change filepath accordingly:
  save(tr_uv05_raw, tr_uv05_hp001, tr_uv05_0120, tr_uv05_011, tr_uv05_515, tr_uv05_15, 
       file = "Z:/non_events/NEvent_2/NEvent_Piton_122009_traces.RData")
  
}else{
  
  # Read in unfiltered data for header:
  uv05_raw <- read.table(file = "Z:/training_data/Train_Piton_012010/UV05_unfiltered.txt", header = TRUE)
  
  # Read in frequency-filtered data to save as trace class:  
  # ======= Change filepath accordingly (repeat for different frequency bands):
  uv05_unfilt <- read.table(file = "Z:/training_data/Train_Piton_012010/UV05_1-5Hz_new.txt", header = TRUE) 
  
  colnames(uv05_unfilt) <- colnames(uv05_raw)
  
  # Convert unix timestamp to date and time. Set to local Reunion time zone = GMT+4.
  
  uv05_unfilt$DateTime <- as.POSIXct(uv05_unfilt$Time, origin="1970-01-01", tz = "GMT")
  
  # Ensure sequence has even length because FFT is very slow on odd length sequences:
  
  if(mod(nrow(uv05_unfilt), 2) == 1){
    
    # Remove first row:
    
    uv05_unfilt <- uv05_unfilt[-1, ]
    
  }
  
  # Replace demo data with Piton data
  #(Note: have not updated network, station, channel, quality, calib, lat/long, elevation, depth, azimut, dip):
  
  tr_uv05 <- tr1
  
  tr_uv05@id <- "Piton"
  tr_uv05@Sensor <- "Unknown sensor"
  tr_uv05@data <- uv05_unfilt$Data
  
  tr_uv05@stats@sampling_rate <- 100 # # 100 Hz (100 readings per second).
  tr_uv05@stats@delta <- 0.01
  tr_uv05@stats@starttime <- uv05_unfilt$DateTime[1]
  tr_uv05@stats@endtime <- uv05_unfilt$DateTime[nrow(uv05_unfilt)]
  tr_uv05@stats@npts <- nrow(uv05_unfilt)
  tr_uv05@stats
  
}

## 3. Compute trace envelopes of frequency-filtered seismic traces.

# ======= Choose the frequency band here.
tr_uv05 <- tr_uv05_011 

piton_ddt <- tr_uv05

na_id <- which(is.na(piton_ddt@data))

temp.time <- proc.time()[3]

if(length(na_id)>0){
  
  nonna_id <- which(!is.na(piton_ddt@data))
  piton_hilbert <- hilbertFFT(piton_ddt@data[nonna_id])
  piton_hilbert_2 <- rep(NA, length(piton_ddt@data))
  piton_hilbert_2[nonna_id] <- piton_hilbert
  sum(is.na(piton_hilbert_2)) == length(na_id) # Check
  piton_hilbert <- piton_hilbert_2
  
}else{
  
  piton_hilbert <- hilbertFFT(piton_ddt@data)
  
}

piton_env_data <- Mod(piton_hilbert) 

piton_env <- tr_uv05
piton_env@id <- "Envelope"
piton_env@data <- piton_env_data 

# Convert to decibel scale:

piton_env_dB <- piton_env
piton_env_dB@data <- 20*log10(piton_env_dB@data)
piton_env_dB@id <- "dB Envelope"

# Check histogram of envelope:

par(mfrow = c(1, 2))
hist(piton_env@data)
hist(piton_env_dB@data)

# Rmd figures: height = 6, width = 8. 
# Visualise for Jan 2010 event (Training set 3): 

if(!miniseed){
  
  s_crisis_start <- as.POSIXct("2010-01-02 07:50:00", tz = "GMT")
  s_swarm_start <- as.POSIXct("2010-01-02 08:20:00", tz = "GMT")
  s_swarm_end <- as.POSIXct("2010-01-02 09:02:00", tz = "GMT")
  eruption_onset_start <- as.POSIXct("2010-01-02 10:20:00", tz = "GMT")
  
  piton_env_dB_2 <- slice(piton_env_dB,  s_crisis_start - 5*60, eruption_onset_start + 180*60)
  
  png(file = "G:/Connected_Extremes/Graphics/train3_envelope.png", width = 2200, height = 2200, res = 300)
  layout(matrix(seq(2)))        # layout a 2x1 matrix
  plot(piton_env_dB, ylab = "Envelope (decibels)", xlab = "Date")
  plot(piton_env_dB_2, ylab = "Envelope (decibels)", xlab = "Time (UTC) on Jan 02")
  abline(v=c(s_crisis_start, s_swarm_start, s_swarm_end,
             eruption_onset_start), col='blue', lwd=2, 
         lty = c(3, 2, 2, 1))
  dev.off()
  
}

# ======= Choose the frequency band here.
save(piton_env_dB, file = "Z:/non_events/NEvent_3/eruption_index_011.RData")

