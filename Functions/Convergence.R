## Plotting and assessing the convergence
### INPUT needed from Bayesian_regression formula
convergencecheck <- function(mcmc_chains, # enter the output from the bayesian_reg function, where you get the sampled values for each chain called xxx$$sampled_values_chains
                             maxlags, # specify a wanted number of lags for the autocorrelation
                             seed # set a seed for reproducibility
){ # with this function you get directly the convergence information regarding trace plots, density plots, autocorrelation and acceptance rates, which are needed to evaluate the convergence of the chains
  set.seed(seed)
  npar <- ncol(mcmc_chains[[1]]) # defining objects which are needed in the function to have it accordingly for the number of parameters chains and iterations
  nchains <- length(mcmc_chains)
  niter <- length(mcmc_chains[[1]][,1]) 
  cols <- brewer.pal(n = 8, name = "Dark2")
  
  # NEW FUNCTION: TRACEPLOTS
  traceplot <- function(mcmc_chains){
    par(mfrow=c(2,3))
    # Plotting the trace plot for all betas - from Intercept (ß0) to number of used predictors
    for(i in 2:npar-1){
      plot(mcmc_chains[[1]][,i],type='l',ylab=substitute(beta[x],list(x=(i-1))),main=substitute(beta[x],list(x=(i-1))),col = cols[1])
      for(h in 2:nchains){
        lines(mcmc_chains[[h]][,i],col=cols[h])} # getting different lines in different colors for each chain
    }
    # Plotting the trace plot for sigma
    plot(mcmc_chains[[1]][,npar],type='l',ylab=expression(sigma^2),main=expression(sigma^2),col=cols[1])
    for(h in 2:nchains){
      lines(mcmc_chains[[h]][,npar],col=cols[h])} # getting different lines in different colors for each chain
  }
  # NEW FUNCTION: DENSITY PLOTS (This is included in the function as I initially planned to include it, but did not have enough space in the report)
  densityplots <- function(mcmc_chains){
    par(mfrow=c(2,3))
    # Plotting the density for all betas - from Intercept (ß0) to number of used predictors
    for(i in 2:npar-1){
      plot(density(mcmc_chains[[1]][,i]),type='l',ylab=substitute(beta[x],list(x=(i-1))),main=substitute(beta[x],list(x=(i-1))),col = cols[1])
      for(h in 2:nchains){
        lines(density(mcmc_chains[[h]][,i]),col=cols[h])} # getting different lines in different colors for each chain
    }
    # Plotting the density for sigma
    plot(density(mcmc_chains[[1]][,npar]),type='l',ylab=expression(sigma^2),main=expression(sigma^2),col=cols[1])
    for(j in 2:nchains){
      lines(density(mcmc_chains[[h]][,npar]),col=cols[h])}
  }
  # NEW FUNCTION: AUTOCORRELATION
  autocorrelation <- function(mcmc_chains, max_lag){
    autocor <- vector("list",nchains) # getting different lines in different colors for each chain
    par(mfrow=c(2,3))
    # obtaining the autocorrelation over all chains and all estimates for a specified number of lags
    for(h in 1:nchains){
      autocor[[h]] <- matrix(0,max_lag,npar)
      
      for(i in 1:max_lag){
        for(j in 1:npar){
          autocor[[h]][i,j] <- cor(mcmc_chains[[h]][1:(niter-(max_lag-1)),j], mcmc_chains[[h]][i:(niter-max_lag+i),j]) # correlation is based on the sampled estimates and the number of specified lags
        }
      }
    }
    # Plotting the autocorrelation of the intercept
    plot(autocor[[1]][,1],type="h",ylab="Intercept",xlab="Lag",col=cols[1])
    for(i in 2:nchains){
      lines(autocor[[i]][,1],type="h",col=cols[i]) # getting different bars in different colors for each chain
    }
    # Plotting the autocorrelation of the predictors
    for(j in 2:(npar-1)){
      plot(autocor[[1]][,j],type="h",ylab=substitute(beta[x],list(x=(j-1))),xlab="Lag",col=cols[1])
      for(i in 2:nchains){
        lines(autocor[[i]][,j],type="h",col=cols[i]) # getting different bars in different colors for each chain
      }
    }
    # Plotting the autocorrelation of the variance
    plot(autocor[[1]][,max(npar)],type="h",ylab=expression(sigma^2),xlab="Lag",col=cols[1])
    for(i in 2:nchains){
      lines(autocor[[i]][,max(npar)],type="h",col=cols[i]) # getting different bars in different colors for each chain
    }
  }
  # Assess the acceptance ratio per chain (was also planned to be included in the report per chain, but not enough space. Only added the pooled acceptance ratio, which can be calculated by the mean of all three chains)
  acceptance_ratio <- vector("list",nchains)
  for(h in 1:nchains){
    acceptance_ratio[[h]] <- matrix(0,npar,1)
    for(k in 1:npar){
      acceptance_ratio[[h]][k] <- round(length(unique(mcmc_chains[[h]][,k]))/(nrow(mcmc_chains[[h]])), digits=3)
    }
  }
  return(list(traceplots=traceplot(mcmc_chains), densityplots=densityplots(mcmc_chains), autocorrelation=autocorrelation(mcmc_chains,maxlags), acceptance_chains= acceptance_ratio))
}