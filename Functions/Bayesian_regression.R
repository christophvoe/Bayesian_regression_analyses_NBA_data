## Bayesian Regression function
### Needed for all other functions as the output of this function is the input of all other functions
bayesian_reg <- function(y, # enter the y variable - labeled before
                         x, # enter the right data set of predictors - all labeled before
                         init, # enter matrix with initial values which fits to the number of predictors in x
                         sigma0, # enter matrix with priors for sigma which fits to the number of predictors in x
                         mu0, # enter matrix with priors for mu which fits to the number of predictors in x
                         burnin, # specify a number of burned samples
                         niter, # specify a number of wanted samples (take number of burn in into account)
                         nchain, # specify a number of chains (take into account that you have to change the initial values matrix for that)
                         seed # set a seed for reproducibility
){ # with this function you can obtain the estimates of a regression procedure using mcmc sampling methods. On top of that you can specify priors and also get an information criteria to compare models
  # even though I don't need the option for one predictor for my analyses I included it as I initially planned to also test one predictor models
  set.seed(seed)  # set a seed
  
  if(is.vector(x)){ # assess the number of predictors in the model. As in case of 1 predictor there is no matrix I included the if else statement
    npred <- ncol(matrix(x))
  } else {
    npred <- ncol(x)
  }
  # create storage for objects including estimates for multiple chains
  samples <- vector("list", nchain) # all sampled estimates (including burn-in)
  samples.burn <- vector("list", nchain) # sampled estimates without burn-in period
  results <- vector("list", nchain) # summary of regression estimates
  
  # create matrix to get pooled estimates over all chains and assigning row and column names to make the output pretty
  results.pooled <- matrix(0,(npred+2),9) 
  if(npred == 1){rownames(results.pooled) <- c("Intercept","Beta 1","Variance")}
  else if(npred == 2){rownames(results.pooled) <- c("Intercept","Beta 1","Beta 2","Variance")}
  else if(npred == 3){rownames(results.pooled) <- c("Intercept","Beta 1","Beta 2","Beta 3","Variance")}
  else {rownames(results.pooled) <- c("Intercept","Beta 1","Beta 2","Beta 3","Beta 4", "Variance")}
  colnames(results.pooled) <- c("Mean","S.E.","MC-error","2.5%","Median","97.5%","Acceptance","Burn-in","Iterations")
  
  for(h in 1:nchain){ # loop over the number of chains
    
    # assign the initial values to two matrices in order to not confuse the variance and the betas
    b <- matrix(init[h,1:(ncol(init)-1)],1,ncol(init)-1) 
    var <- matrix(init[h,ncol(init)],1,1)
    
    # create one matrix per chain
    samples[[h]] <- matrix(0,niter,(npred+2))
    samples[[h]][1,] <- c(b,var)
    results[[h]] <- matrix(0,(npred+2),9)
    
    # for loop - sample as often as iterations are specified in the function
    for(i in 2:niter){
      # Gibbs step to sample the intercept
      V <- 1/(length(y)/var + 1/sigma0[1])
      if(npred == 1){
        M <- (sum(y - x*b[,-1])/var + mu0[1]/sigma0[1])*V
      } else {
        M <- (sum(y - x%*%b[,-1])/var + mu0[1]/sigma0[1])*V # for more then 1 predictor using matrix algebra as we need all x's with the corresponding beta, except the intercept <- (b[,-1])
      }
      b[,1] <- rnorm(1,M,sqrt(V)) # sampling the value for beta 0 and save it in matrix b in first column
      
      # Metropolis Hastings step - Sampling beta 1 with the Metropolis Hastings step
      b.curr <- b[,2]
      b.prop <- rnorm(1,b.curr,0.005) # parameter which can be used to obtain changes regarding the acceptance ratio and autocorrelation - I already played around with it and I am most satisfied with this value
      u <- runif(1,0,1) # sample random value between 0 and 1, which will be compared to my proposed estimate
      
      if(npred == 1){ # again adapting the way to get my likelihood based on the number of used predictors
        V <- 1/(sum(x^2)/var)
        M <- sum(x*(y - b[,1]))/var*V
      } else if(npred == 2){
        V <- 1/(sum(x[,1]^2)/var)
        M <- sum(x[,1]*(y - b[,1] - x[,-1]*b[,-c(1:2)]))/var*V
      } else {
        V <- 1/(sum(x[,1]^2)/var)
        M <- sum(x[,1]*(y - b[,1] - x[,-1]%*%b[,-c(1:2)]))/var*V
      }
      
      likelihood.curr <- dnorm(b.curr,M,sqrt(V),log=TRUE) # obtaining the likelihood of the current and the proposed value based on the calculations above
      likelihood.prop <- dnorm(b.prop,M,sqrt(V),log=TRUE)
      
      prior.curr <- dnorm(b.curr,mu0[2],sigma0[2],log=TRUE) # obtaining the likelihoods of the current and the proposed priors, which are needed for the calculations
      prior.prop <- dnorm(b.prop,mu0[2],sigma0[2],log=TRUE)
      
      target.curr <- likelihood.curr + prior.curr # Again we need for the current and proposed option # as we are using the log transformed values the signs change (+ instead of *) - I am doing this as my data might be too skewed and this was proposed as one way to get rid of any estimation problems
      target.prop <- likelihood.prop + prior.prop
      
      if(exp(target.prop - target.curr) > u){ # this statement evaluates if our proposed value is higher then the sampled value in u (between 0 and 1) - if it is higher we take our proposed estimate as the newest estimate otherwise we take the old value again as the estimate
        b[,2] <- b.prop
      } else {
        b[,2] <- b.curr
      }
      
      # Sampling all further betas again with Gibbs sampling - However, right now the function only works with three predictors as I specified some things above and below in a certain way. Could be altered for the use of infinetly many predictors
      if(npred > 1){
        for(j in 2:npred){ # this for loop runs from 2 to the specified number of predictors
          V <- 1/(sum(x[,j]^2)/var + 1/sigma0[(j+1)])
          if(npred == 2){
            M <- (sum(x[,j]*(y - b[,1] - x[,-j]*b[,-c(1,j+1)]))/var + mu0[(j+1)]/sigma0[(j+1)])*V
          } else {
            M <- (sum(x[,j]*(y - b[,1] - x[,-j]%*%b[,-c(1,j+1)]))/var +
                    mu0[(j+1)]/sigma0[(j+1)])*V
          }
          b[,(j+1)] <- rnorm(1,M,sqrt(V)) # rest is the same compared to other Gibbs sampling steps
        }
      }
      
      # Sampling the variance
      A <- length(y)/2 + alpha0
      if(npred == 1){
        B <- sum((y - (b[,1] + x*b[,-1]))^2)/2 + beta0
      } else {
        B <- sum((y - (b[,1] + x%*%b[,-1]))^2)/2 + beta0
      }
      var <- 1/rgamma(1,shape=A,rate=B) # here we use the gamma function and taking the inverse of it. This is different to the other sampled values, but needed to get the estimates we want!
      
      samples[[h]][i,] <- c(b,var) # combine the two data sets of sampled values and save them in one data set.
    } 
    samples.burn[[h]] <- matrix(0,(niter-burnin),(npred+2)) # making a matrix where we can save only the values we are interested - so we  discard burn-in period
    samples.burn[[h]] <- samples[[h]][(burnin+1):niter,]
    samples.all <- do.call(rbind, samples.burn)  # combine the samples of each chain in one data set - In this way we can obtain the pooled estimates
    
    for(r in 1:(npred+2)){ # calculate the estimates, which appear in the output based on the pooled estimates
      results.pooled[r,1] <- round(mean(samples.all[,r]), digits=3)
      results.pooled[r,2] <- round(sd(samples.all[,r]), digits=3)
      results.pooled[r,3] <- round(sd(samples.all[,r])/sqrt(niter*nchain-burnin*nchain), digits=5)
      results.pooled[r,4] <- round(quantile(samples.all[,r],0.025), digits=3)
      results.pooled[r,5] <- round(median(samples.all[,r]), digits=3)
      results.pooled[r,6] <- round(quantile(samples.all[,r],0.975), digits=3)
      results.pooled[r,7] <- round((length(unique(samples.all[,r]))/(niter-burnin))/nchain,digits=3)
      results.pooled[r,8] <- burnin*nchain
      results.pooled[r,9] <- niter*nchain-burnin*nchain
    }
    # Integrating the DIC step in the function as we otherwise have to specify same information again in a new function and as we try to model a regression - an information criteria should be included in the function itself as well
    twologlike <- rep(NA, nrow(samples.all)) # I am separating the code for the different amount of predictors again due to simplicity 
    if(npred==1){  # DIC for one predictor
      for(i in 1:nrow(samples.all)){ # calculating the likelihood of each iteration for dbar
        twologlike[i] <- -2*sum(dnorm(y, mean = samples.all[i,1] + samples.all[i,2]*x, sd=sqrt(samples.all[i,3]), log = TRUE))}
      dbar <- mean(twologlike) # Take the mean over the separate likelihoods
      dhat <- -2*sum(dnorm(y, mean = results.pooled[1,1] + results.pooled[2,1]*x, sd=sqrt(results.pooled[3,1]), log = TRUE)) # Calculate Dhat (calculate the -2 log likelihood with the posterior means of the parameters)
    }
    else if(npred==2){ # DIC for two predictors
      for(i in 1:nrow(samples.all)){
        twologlike[i] <- -2*sum(dnorm(y, mean = samples.all[i,1] + samples.all[i,2]*x[,1] + samples.all[i,3]*x[,2], sd=sqrt(samples.all[i,4]), log = TRUE))}
      dbar <- mean(twologlike)
      dhat <- -2*sum(dnorm(y, mean = results.pooled[1,1] + results.pooled[2,1]*x[,1] + results.pooled[3,1]*x[,2], sd = sqrt(results.pooled[4,1]), log = TRUE))
    }
    else{   # DIC for four predictors
      for(i in 1:nrow(samples.all)){
        twologlike[i] <- -2*sum(dnorm(y, mean = samples.all[i,1] + samples.all[i,2]*x[,1] + samples.all[i,3]*x[,2] + samples.all[i,4]*x[,3] + samples.all[i,5]*x[,4], sd=sqrt(samples.all[i,6]), log = TRUE))}
      dbar <- mean(twologlike)
      dhat <- -2*sum(dnorm(y, mean = results.pooled[1,1] + results.pooled[2,1]*x[,1] + results.pooled[3,1]*x[,2] +  results.pooled[4,1]*x[,3]  +  results.pooled[5,1]*x[,4], sd = sqrt(results.pooled[6,1]), log = TRUE))
    }
    dic <- dhat + 2*(dbar - dhat) # object that can be used to compare different models
  }
  print(results.pooled) # give pooled estimates as default output - Rest can be obtained by saving an object where you use the function and ask for specific output
  return(list(estimates=results.pooled, estimates_chains=results, sampled_values=samples.all, sampled_values_chains=samples.burn, burnin=samples[1:burnin], DIC = dic))
}