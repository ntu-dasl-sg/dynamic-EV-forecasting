##################################################################
## This R script computes the area under the curve and other    ##
## evaluation metrics to evaluate the value of using extreme    ##
## value theory to choose the threshold for the dynamic extreme ##
## value model based on the January 2010 eruption data          ##
## (Training Set 3).                                            ##
################################################################## 

library(evd)
library(ggplot2)
library(MASS)
library(gridExtra)
library(car)
library(xtable)
library(generalhoslem) # Hosmer-Lemeshow Test
library(cutpointr) # ROC and optimal cutoff points.
library(SpecsVerification) # ReliabilityDiagram.
library(tidyr) # For reshaping dataframes.

setwd("D:/Documents/Imperial_NTU_collaboration/Seismic data")

source('D:/Documents/Imperial_NTU_collaboration/Seismic data/GPD_regression.r')

freq_bands <- c("011", "15", "515", "0120", "hp001")

## 1. Set thresholds selected by extreme value theory for the respective frequency-filtered envelopes.
ev_thres <- c(69, 85, 84, 88, 89)
zero_shape_v <- c(TRUE, FALSE, FALSE, TRUE, TRUE)

## 2. First run with thresholds set at chosen values
# Create starting dataframe:
# accuracy_df <- data.frame("freq_band" = NA, "ev_thres" = NA, "no_exceed" = NA, 
#                           "sensitivity" = NA, "specificity" = NA, 
#                           "AUC" = NA, "brier" = NA, "optimal_c" = NA, 
#                           "reliab_x1" = NA, "reliab_x2" = NA, "reliab_x3" = NA, 
#                           "reliab_x4" = NA, "reliab_x5" = NA, "reliab_x6" = NA,
#                           "reliab_x7" = NA, "reliab_x8" = NA, "reliab_x9" = NA, 
#                           "reliab_x10" = NA,
#                           "reliab_y1" = NA, "reliab_y2" = NA, "reliab_y3" = NA, 
#                           "reliab_y4" = NA, "reliab_y5" = NA, "reliab_y6" = NA,
#                           "reliab_y7" = NA, "reliab_y8" = NA, "reliab_y9" = NA, 
#                           "reliab_y10" = NA)
# Run for loop below with run_j = 1 only. 

## 3. Subsequent runs for run_j = 0.9, 0.7, 0.5.

# Read in first run results.
accuracy_df <- read.csv("accuracy_df.csv")

temp.time <- proc.time()[3]

for (run_j in seq(0.9, 0.5, by = -0.1)){ # thres_pc
  
  for (run_i in 1:5){ # freq_bands
    
    ## Run settings:
    
    freq_i <- run_i 
    thres_pc <- run_j 
    
    zero_shape <- zero_shape_v[freq_i]
    
    # Read in data:
    
    load(file = paste("Piton/index_trace_env_", freq_bands[freq_i], "/lag1h_model_df_", freq_bands[freq_i], ".RData", sep = "")) 
    
    model_df <- lag1h_model_df
    
    eruption_index <- model_df$piton_env_dB
    
    u_chosen <- thres_pc*ev_thres[freq_i]
    
    # Times of exceedance
    exceed_id <- which(model_df$piton_env_dB >= u_chosen)
    exceed_times <- model_df$DateTime[exceed_id]
    
    model_df$exceed <- 0
    model_df$exceed[exceed_id] <- 1
    
    ## Covariate transformations
    
    cov_col <- which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB", "exceed")))
    length(cov_col) # 280.
    cov_names <- colnames(model_df)[cov_col]
    
    # Box-Cox transformation for more normal distribution shapes:
    
    lambda_vector <- rep(NA, length(cov_col))
    buffer_vector <- rep(NA, length(cov_col))
    
    test_col <- c(1:length(cov_col))
    
    
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
      sd_vector[i] <-temp_sd
      
      stand_cov <- (covariate - temp_mean)/temp_sd
      
      model_df[, cov_col[i]] <- stand_cov
      
    }
    
    trans_df$mean <- mean_vector
    trans_df$sd <- sd_vector
    
    nan_col <- which(lapply(model_df, FUN = function(x){sum(is.nan(x))}) > 0)
    
    model_df <- model_df[, !(colnames(model_df) %in% names(nan_col))]
    
    ## Logistic regression for threshold exceedances
    
    te_df <- model_df[, which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB")))]
    
    # Order covariates in terms of how much they inform univariate models:
    
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
    
    te_stepAIC <- stepAIC(te_logistic, direction = "both")
    
    te_model <- te_stepAIC
    te_summary <- summary(te_model)
    
    te_lp_predict <- predict(te_model, type = "link", se.fit = TRUE)
    
    model_df$te_lp <- te_lp_predict$fit
    model_df$te_lp_se <- te_lp_predict$se.fit
    model_df$te_lp_upper <- te_lp_predict$fit + qnorm(0.975)*te_lp_predict$se.fit
    model_df$te_lp_lower <- te_lp_predict$fit - qnorm(0.975)*te_lp_predict$se.fit
    
    model_df$te_prob <- 1 / (1 + exp(-model_df$te_lp)) # Inverse logit.
    model_df$te_prob_upper <- 1 / (1 + exp(-model_df$te_lp_upper))
    model_df$te_prob_lower <- 1 / (1 + exp(-model_df$te_lp_lower))
    
    # Evaluating forecast accuracy
    
    # Receiver Operating Characteristic: 
    
    cp <- cutpointr(model_df, te_prob, exceed, 
                    method = maximize_metric, metric = youden)
    
    reliability_diag <- ReliabilityDiagram(
      model_df$te_prob,
      model_df$exceed,
      bins = 10,
      nboot = 500,
      plot = FALSE,
      plot.refin = FALSE,
      cons.probs = 0.95,
      attributes = FALSE,
      handle.na = c("na.fail", "use.pairwise.complete")
    )
    
    # Brier score:
    brier_score <- mean((model_df$te_prob - model_df$exceed)^2)
    
    # Save accuracy results:
    
    accuracy_v <-c(freq_bands[freq_i], u_chosen, length(exceed_times), 
                   cp$sensitivity, cp$specificity, 
                   cp$AUC, brier_score, cp$optimal_cutpoint, 
                   reliability_diag$p.avgs,
                   reliability_diag$cond.probs)
    
    accuracy_df <- rbind(accuracy_df, accuracy_v)
    
    write.csv(accuracy_df, "accuracy_df.csv", row.names = FALSE)
    
    ## Dynamic GP model for threshold excesses
    
    test_model_1 <- fpot(model_df$piton_env_dB, threshold = u_chosen, model = "gpd", cmax = FALSE)
    test_model_1
    fixed_shape <- test_model_1$estimate['shape'] 
    
    model_df$excess <- model_df$piton_env_dB - u_chosen 
    
    excess_df <- model_df[, which(!(colnames(model_df) %in% c("Time", "Data", "DateTime", "piton_env", "piton_env_dB", "piton_env_dB_runmed", "te_prob", "exceed")))]
    
    excess_df <- excess_df[excess_df$excess>=0, ]
    
    # Order covariates in terms of how much they inform univariate models:
    
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
    
    
    if(zero_shape){
      
      excess_expreg <- glm(excess ~ ., family = Gamma(link = "log"), data = excess_df)
      # 10 not defined because of singularities => probably should have check for multicollinearity first.
      
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
    
    png(paste("Graphics/gpdfit_", freq_bands[freq_i], "_", thres_pc, ".png", sep = ""), width = 2000, height = 2000, res = 300)
    
    par(mfrow = c(2, 2))
    plot(cgp_model) # Q-Q plot uses model based estimate of the quantile function and the data used in the fitted model, sorted into ascending order.
    
    dev.off()
    
    
  }
  
}

time.taken <- proc.time()[3] - temp.time # 14 minutes for 1 run.

## 4. Plot changes in AUC with threshold choice.

accuracy_df <- read.csv("Y:/accuracy_df.csv")

freq_bands_2 <- unique(accuracy_df$freq_band) # csv cannot read characters starting with "0".

for (i in freq_bands_2){
  
  temp_df <- accuracy_df[accuracy_df$freq_band == i, ]
  temp_df <- temp_df[order(temp_df$ev_thres, decreasing = TRUE), ]
  accuracy_df[accuracy_df$freq_band == i, ] <- temp_df
  
}

accuracy_df$thres_pc <- rep(c(1, 0.9, 0.8, 0.7, 0.6, 0.5), each = 5)

accuracy_df$freq_band <- factor(accuracy_df$freq_band, levels = c("11", "15", "515", "120", "hp001"))
levels(accuracy_df$freq_band) <- c("0.1-1 Hz", "1-5 Hz", "5-15 Hz", "0.1-20 Hz", "High pass 0.01Hz")

accuracy_df_trunc <- accuracy_df[accuracy_df$freq_band %in% c("1-5 Hz", "0.1-20 Hz", "High pass 0.01Hz"), ]

png("D:/Documents/Imperial_NTU_collaboration/Seismic data/Graphics/AUC_thres_pc_all_freq_trunc.png", width = 1500, height = 1250, res = 300)

ggplot(data = accuracy_df_trunc) + geom_line(aes(x = thres_pc, y = AUC, linetype = freq_band)) + geom_point(aes(x = thres_pc, y = AUC, shape = freq_band)) + theme_classic() + labs(x = "Proportion of EV threshold", shape = "Frequency band", linetype = "Frequency band") + theme(legend.position="bottom")

dev.off()

