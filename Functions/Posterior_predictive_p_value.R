## Posterior predictive p-value
### Two self-assessed metrics are used
### INPUT needed from Bayesian_regression formula

ppp <- function(sampled_values, # enter the output from the bayesian_reg function, where you get the sampled values for each chain called xxx$$sampled_values
                y, # enter the y variable - labeled before 
                x, # enter the right data set of predictors - all labeled before
                nsims, # enter the number of wanted simulated data sets
                seed # set a seed for reproducibility
){ # with this function I calculate the residuals of simulated data sets and use them and the fitted values to execute two ppp checks. As we simulate normal data our observed data set should be somewhere around 0.5 for both ppp checks to have a fitting model.
  set.seed(seed)
  if(is.vector(x)){
    npred <- ncol(matrix(x))
  } else {
    npred <- ncol(x) # assessing the number of used predictors in the model
  }
  # storage for all the different needed objects
  simdata <- array(data = NA, dim = c(nsims, length(y)))
  fitsim <- array(data = NA, dim = c(nsims, length(y)))
  corsimfit <- rep(0, nsims)
  obsresiduals <- array(data = NA, dim = c(nsims, length(y)))
  corobsfit <- rep(0, nsims)
  pvals.hom <- rep(NA, nsims)
  # storage for second test
  distsimfit <- rep(0, nsims)
  distobsfit <- rep(0, nsims)
  pvals.lin <- rep(NA, nsims)
  
  # simulating the residuals and the fitted values according to the number of predictors
  if(npred==2){ # does the same for two predictors
    for(i in 1:nsims){
      simdata[i,] <- rnorm(length(y), mean = sampled_values[i,1] + x[,1]*sampled_values[i,2] + x[,2]*sampled_values[i,3], sd = sqrt(sampled_values[i,4]))
    }
    for(i in 1:nsims){
      fitsim[i,] <- sampled_values[i,1] + x[,1]*sampled_values[i,2] + x[,2]*sampled_values[i,3]
    }
    simresiduals <- simdata -fitsim
  }
  else{ # does the same for four predictors
    for(i in 1:nsims){
      simdata[i,] <- rnorm(length(y), mean = sampled_values[i,1] + x[,1]*sampled_values[i,2] + x[,2]*sampled_values[i,3] + x[,3]*sampled_values[i,4]+ x[,4]*sampled_values[i,5], sd = sqrt(sampled_values[i,6]))
    }
    for(i in 1:nsims){
      fitsim[i,] <- sampled_values[i,1] + x[,1]*sampled_values[i,2] + x[,2]*sampled_values[i,3] + x[,3]*sampled_values[i,4] + x[,4]*sampled_values[i,5]
    }
    simresiduals <- simdata -fitsim
  }
  # first test statistic regarding homoscedasticity - correlation between fitted and residual values
  for(i in 1:nsims){
    corsimfit[i] <- cor(simresiduals[i,], fitsim[i,])
  }
  for(i in 1:nsims){
    obsresiduals[i,] <- y - fitsim[i,]
  }
  for(i in 1:nsims){
    corobsfit[i] <- cor(obsresiduals[i,], fitsim[i,])
  }
  # second test regarding linearity - mean value of distance between residual and fitted value
  for(i in 1:nsims){
    distsimfit[i] <- mean(simresiduals[i,1:(length(y)*0.33)]) + mean(simresiduals[i,((length(y)*0.33)+1):(length(y)*0.66)]) + mean(simresiduals[i,((length(y)*0.66)+1):length(y)]) # I split the distribution of residuals in three equal parts and take the mean of the residuals in each part. Next, I take the sum of the three parts to get one single value for each simulated data set
  }
  
  for(i in 1:nsims){
    distobsfit[i] <-  mean(obsresiduals[i,1:(length(y)*0.33)]) + mean(obsresiduals[i,((length(y)*0.33)+1):(length(y)*0.66)]) + mean(obsresiduals[i,((length(y)*0.66)+1):length(y)]) # here I do the same for the observed residuals
  }
  for(i in 1:nsims){
    pvals.hom[i] <- ifelse(corsimfit[i] > corobsfit[i], 1, 0) # evaluates the proportion of times the simulated data set has a higher correlation then the observed dataset
  }
  pvals.hom <- round(mean(pvals.hom),3) # gets us the PPP
  
  for(i in 1:nsims){
    pvals.lin[i] <- ifelse(distsimfit[i] > distobsfit[i], 1, 0) # evaluates the proportion of times the simulated data set has a higher sum then the observed dataset
  }
  pvals.lin <- round(mean(pvals.lin),3) # gets us the PPP
  
  return(list(Homoscedactisty=pvals.hom, Linearity=pvals.lin))
}