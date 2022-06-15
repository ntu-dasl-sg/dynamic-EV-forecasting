##################################################################
## These R functions compute the likelihood of the Generalised  ##
## Pareto Distribution (GPD) regression as in the Bee et al.    ##
## (2019) paper. They feed into the GPD_stepAIC functions.       ##
################################################################## 

# The GPD_stepAIC function defines X (covariate matrix) and z_t (excess).

GPD_nll <- function(beta, fixed_shape, X, z_t){
  
  v_t <- exp(X%*%beta) # dimension: 544 x 1.
  
  T2 <- 1 + (fixed_shape/v_t)*z_t
  if(sum(T2<=0)>0){nll <- 99999}else{
    power_val <- -(1/fixed_shape)-1
    nll <- -sum(-log(v_t) + power_val*log(T2))
  }
  
  return(nll)

}

# To re-estimate shape again:
GPD_nll_2 <- function(beta, X, z_t){
  
  fixed_shape <- beta[1]
  v_t <- exp(X%*%beta[-1]) # dimension: 544 x 1.
  
  T2 <- 1 + (fixed_shape/v_t)*z_t
  if(sum(T2<=0)>0){nll <- 99999}else{
    power_val <- -(1/fixed_shape)-1
    nll <- -sum(-log(v_t) + power_val*log(T2))
  }
  
  return(nll)
  
}

##################################################################
## This R function conducts a step-wise variable selection for  ##
## the Generalised Pareto Distribution (GPD) regression based   ##
## minimising the Akaike Information Criterion (AIC). It        ##
## follows the default direction of forwards then backwards     ## 
## selection.                                                   ##
################################################################## 

GPD_stepAIC <- function(excess_df, fixed_shape, max_it = 100, fix = TRUE){

cov_col <- which(!(colnames(excess_df) %in% c('exceed', 'excess')))
X_full <- cbind(matrix(1, nrow = nrow(excess_df), ncol = 1), as.matrix(excess_df[, cov_col]))
z_t <- excess_df$excess

iter_round <- 1
current_cov <- NULL
available_cov <- cov_col

AIC_store <- NULL

while(iter_round < max_it){
  
  round_AIC <- rep(NA, length(available_cov))
  
  for(i in 1:length(available_cov)){
    
    test_cov <- c(current_cov, available_cov[i])
    X <- X_full[, c(1, 1+test_cov)]
    start_beta <- c(5, rep(0, ncol(X)-1)) # Intercept initial value about log(150).

    GPD_MLE <- optim(start_beta, GPD_nll, X = X, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE)
    
    round_AIC[i] <- 2*GPD_MLE$value + 2*length(GPD_MLE$par)
    
  }
  
  chosen_id <- order(round_AIC)[1] # Index in available_cov terms: 1-length(available_cov).
  
   new_AIC <- round_AIC[chosen_id]
  
  if(new_AIC<min(AIC_store)){
  
    AIC_store <- c(AIC_store, round_AIC[chosen_id])  
    
    chosen_cov <- available_cov[chosen_id] # Index in cov_col terms.
    
    current_cov <- c(current_cov, chosen_cov) # Index in cov_col terms.
    
    available_cov <- available_cov[-chosen_id] # Index in cov_col terms.
    
    iter_round <- iter_round + 1  
    
  }else{
    
    
    print("No further decreases in AIC. Stop stepAIC.")
    break  
    
  }
  
  # Consider backwards selection from iter_round = 4 onwards. (at least 3 current covariates)
  
  if(iter_round > 3){
    
    back_AIC <- rep(NA, length(current_cov))
    
    for(i in 1:length(current_cov)){
      
      test_cov <- current_cov[-i]
      X <- X_full[, c(1, 1+test_cov)]
      start_beta <- c(5, rep(0, ncol(X)-1)) # Intercept initial value about log(150).
 
      GPD_MLE <- optim(start_beta, GPD_nll, X = X, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE) 
      
      back_AIC[i] <- 2*GPD_MLE$value + 2*length(GPD_MLE$par)
      
      
    }
    
    if(min(back_AIC) < min(AIC_store)){
      
      return_cov <- which(back_AIC == min(back_AIC))
      current_cov <- current_cov[-return_cov]
      available_cov <- c(available_cov, current_cov[return_cov])
    }
    
  }
  
  print(paste("Iteration", iter_round, "finished.", sep = " "))
  
}

return(list("current_cov" = current_cov, "AIC_store" = AIC_store))

}

##################################################################
## This R function conducts a step-wise variable selection for  ##
## the Generalised Pareto Distribution (GPD) regression based   ##
## minimising the Akaike Information Criterion (AIC). It        ##
## follows the stepAIC default direction of backward then       ## 
## forward selection.                                           ##
################################################################## 

GPD_stepAIC2 <- function(excess_df, fixed_shape, max_it = 100, fix = TRUE){ 
  
  cov_col <- which(!(colnames(excess_df) %in% c('exceed', 'excess')))
  X_full <- cbind(matrix(1, nrow = nrow(excess_df), ncol = 1), as.matrix(excess_df[, cov_col]))
  z_t <- excess_df$excess
  
  iter_round <- 1
  current_cov <- cov_col
  available_cov <- NULL
  
  # Full model
  start_beta <- c(5, rep(0, ncol(X_full)-1))
  GPD_MLE <- optim(start_beta, GPD_nll, X = X_full, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE)
  
  AIC_store <- 2*GPD_MLE$value + 2*length(GPD_MLE$par)
  
  
  while(iter_round < max_it){
    
    round_AIC <- rep(NA, length(current_cov))
    
    for(i in 1:length(current_cov)){
      
      test_cov <- current_cov[-i]
      X <- X_full[, c(1, 1+test_cov)]
      start_beta <- c(5, rep(0, ncol(X)-1)) # Intercept initial value about log(150).

      GPD_MLE <- optim(start_beta, GPD_nll, X = X, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE) 
      
      round_AIC[i] <- 2*GPD_MLE$value + 2*length(GPD_MLE$par)
      
    }
    
    chosen_id <- order(round_AIC)[1] # Index in current_cov terms: 1-length(current_cov).
    
    new_AIC <- round_AIC[chosen_id]
    
    if(new_AIC<min(AIC_store)){
      
      AIC_store <- c(AIC_store, round_AIC[chosen_id])  

      available_cov <- c(available_cov, current_cov[chosen_id]) # Index in cov_col terms.

      current_cov <- current_cov[-chosen_id] # Index in cov_col terms.
            
      iter_round <- iter_round + 1  
      
    }else{
      
      
      print("No further decreases in AIC. Stop stepAIC.")
      break  
      
    }
    
    # Consider backwards selection from iter_round = 4 onwards. (at least 3 current covariates)
    
    if(iter_round > 3){
      
      back_AIC <- rep(NA, length(available_cov))
      
      for(i in 1:length(available_cov)){
        
        test_cov <- c(current_cov, available_cov[i])
        X <- X_full[, c(1, 1+test_cov)]
        start_beta <- c(5, rep(0, ncol(X)-1)) # Intercept initial value about log(150).

        GPD_MLE <- optim(start_beta, GPD_nll, X = X, z_t = z_t, fixed_shape = fixed_shape, method = "Nelder-Mead", hessian = FALSE)
        
        back_AIC[i] <- 2*GPD_MLE$value + 2*length(GPD_MLE$par)
        
      }
      
      if(min(back_AIC) < min(AIC_store)){
        
        add_cov <- which(back_AIC == min(back_AIC))
        current_cov <- c(current_cov, available_cov[add_cov])
        available_cov <- available_cov[-add_cov]
      }
      
    }
    
    print(paste("Iteration", iter_round, "finished.", sep = " "))
    
  }
  
  return(list("current_cov" = current_cov, "AIC_store" = AIC_store))
  
}


# fisher_info<-solve(GPD_MLE$hessian)
# prop_sigma<-sqrt(diag(fisher_info))
# upper<-GPD_MLE$par+qnorm(0.975)*prop_sigma
# lower<-GPD_MLE$par-qnorm(0.975)*prop_sigma
# interval<-data.frame(value=GPD_MLE$par, lower=lower, upper=upper, se = prop_sigma)


# # Another example
# gdata <- excess_df[, 1:2]
# Threshold <- 0; xi <- fixed_shape
# gdata <- transform(gdata, y2 = excess_df$excess, shape = xi)
# fit <- vglm(y2 ~ ., gpd(Threshold), data = gdata, trace = TRUE)
# coef(fit, matrix = TRUE)
