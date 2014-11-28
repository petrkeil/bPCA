
sim.bPCA <- function(data, 
                    #priors specification
                    covmat.prior,
                    covmat.prior.DF,
                    mu.prior,
                    mu.prior.cov,
                    #MCMC specication
                    n.chains,
                    n.iter,
                    n.burnin){
  
  # requirements
  require(R2jags)
  require(MASS)
  require(Matrix)
  require(coda)
  
  # dataset dimensions
  N = nrow(data)
  V = ncol(data)
  
  # defaults

  if(missing(covmat.prior)) covmat.prior=as.matrix(Diagonal(V,1/1000))
  if(missing(covmat.prior.DF)) covmat.prior.DF=V
  if(missing(mu.prior)) mu.prior=rep(0,V)
  if(missing(mu.prior.cov)) mu.prior.cov=as.matrix(Diagonal(V,1000))
  if(missing(n.chains)) n.chains=3
  if(missing(n.iter)) n.iter=5000
  if(missing(n.burnin)) n.burnin=4500
    
  # makes precisions from covariances
  mu.prior.prec=ginv(mu.prior.cov)
  
  # puts data into list
  listdata = list(Y=as.matrix(data), 
                  N=N,
				          V=V, 
                  covmat.prior=covmat.prior, 
                  mu.prior=mu.prior, 
                  covmat.prior.DF=covmat.prior.DF, 
                  mu.prior.prec=mu.prior.prec)
  
  # defines the model in JAGS language
  cat("
  model
  {
    # priors on the vector of multinormal means
    mu[1:V] ~ dmnorm(mu.prior[], mu.prior.prec[,])
    # priors on the covariance matrix
    prec[1:V,1:V] ~ dwish(covmat.prior[,], covmat.prior.DF)
    
    # makes covariance from precision
    cov[1:V,1:V] <- inverse(prec[,])

    # likelihood
    for (i in 1:N)
    {
      Y[i,1:V] ~ dmnorm(mu[], prec[,])
    }
  }
  ", file="PCA.bugs")
  
  # jags model to estimate covariance matrix distribution
  pcabay <- jags(data=listdata,
                 model.file="PCA.bugs",
                 parameters.to.save=c("cov", "mu"),
                 n.chains=n.chains,
                 n.iter=n.iter,
                 n.burnin=n.burnin, 
                 DIC=FALSE)
  
  return(pcabay)
}

# ------------------------------------------------------------------------------
# CONVERGENCE DIAGNOSTICS

convergenceplots.bPCA <- function(bPCA.fitted)  
{
   #convergence diagnostics plots of the means mu
   plot(as.mcmc(bPCA.fitted))
} 

# ------------------------------------------------------------------------------ 
# BOXPLOTS OF THE "STABILITY" OF THE EIGENVALUES AND EXPLAINED VARIANCES

eigenvalplots.bPCA <- function(bPCA.fitted, data)  
{
    V <- ncol(data)
    sims <- bPCA.fitted$BUGSoutput$sims.matrix
    # extracting only the covariances
    sims <- sims[,1:(V*V)]
    # empty matrix for results
    eigen.chains <- matrix(nrow=nrow(sims), ncol=V)
    
    # calculate eigenvalues for each covariance matrix in the chain
    for(i in 1:nrow(sims))
    {
      covm <- matrix(sims[i,], V, V)
      eigen.chains[i,] <- eigen(covm)$values
    }
    # percents of explained variability
    exp.vars <- eigen.chains/rowSums(eigen.chains) * 100
    # posteriors of eigenvalues as boxplots
    par(mfrow=c(1,2))
    boxplot(eigen.chains, ylab="Eigenvalue", xlab="PCA axis", 
            col="grey", outline=FALSE)
    boxplot(exp.vars, ylab="Explained variability [% of total]", xlab="PCA axis", 
            col="grey", outline=FALSE, ylim=c(0,100))
    
    results <- list(Eigenvalues = summary(eigen.chains),
                    Exp.var = summary(exp.vars))
    return(results)
} 

# ------------------------------------------------------------------------------ 

plot.classicPCA <- function(data, axes.to.plot=1:2, scale=1, xlim, ylim)
{
  eig=eigen(cov(data))
  loadings <- t(t(eig$vectors))*(eig$values^scale)
  row.names(loadings) <- names(data)
  centered <- scale(data, scale=FALSE)
  scores <-  centered  %*% eig$vectors

  scores <- scores[,axes.to.plot]
  loadings <- loadings[,axes.to.plot]
  
  biplot(scores,loadings, main="Classic PCA",
         xlab=paste("Comp.", axes.to.plot[1]),
         ylab=paste("Comp.", axes.to.plot[2]),
         ylim=ylim, xlim=xlim)
  abline(h=0, lty=2, col="grey"); abline(v=0, lty=2, col="grey")
}

 
# ------------------------------------------------------------------------------
# SIMPLE BIPLOTS OF THE BAYESIAN PCA

biplots.bPCA <- function(bPCA.fitted, data, axes.to.plot=1:2, scale=1)
{
  V = length(data[1,])
  summ.stats <- c("2.5%", "50%", "97.5%")
  par(mfrow=c(1, length(summ.stats)))
  
  for(summ.stat in summ.stats)
  {
      covm = matrix(bPCA.fitted$BUGSoutput$summary[1:(V^2), summ.stat], V, V)
      mu = bPCA.fitted$BUGSoutput$summary[((V^2)+1):((V^2)+V), summ.stat]
      eig=eigen(covm)
      loadings <- t(t(eig$vectors)*(eig$values^scale))
      row.names(loadings) <- names(data)
      centered <- scale(data, center=mu, scale=FALSE)
      scores <- centered  %*% eig$vectors
      biplot(x=scores[,axes.to.plot], 
             y=loadings[,axes.to.plot], 
             main=paste(summ.stat, "of Bayesian PCA"), 
             xlab=paste("Comp.", axes.to.plot[1]),
             ylab=paste("Comp.", axes.to.plot[2]))
      abline(h=0, lty=2, col="grey"); abline(v=0, lty=2, col="grey")
  }
}

# ------------------------------------------------------------------------------
# FAST LOOP-AVOIDING FUNCTION THAT EXTRACTS THE CHAINS OF THE LOADINGS

get.loadings.chain.bPCA <- function(bPCA.fitted, data)
{
  sims <- bPCA.fitted$BUGSoutput$sims.matrix
  V <- ncol(data)
  sims.cov <- sims[,1:(V^2)]
  # split by rows (for the lapply function below)
  sims.cov <- split(sims.cov, seq(nrow(sims.cov)))
  names(sims.cov) <- NULL
  
  # the function that will be used by lapply
  load.extract <- function(cov, V, data)
  {
    covm = matrix(cov, V, V)
    loadings <- eigen(covm)$vectors
    row.names(loadings) <- names(data)
    colnames(loadings) <- paste("Comp.", 1:V, sep="")
    return(loadings)
  }
  
  # the loop-avoiding lapply that applies the load.extract() to
  # each element of the sims.cov list
  loadings.chain <- lapply(X=sims.cov, FUN=load.extract, V=V, data=data)
  return(loadings.chain)
}

# ------------------------------------------------------------------------------
# FAST FUNCTION THAT SUMMARIZES THE CHAINS FOR THE LOADINGS BY QUANTILES
# AND BY LATTICE HISTOGRAMS

# plots and summarizes loadings
summary.loadings.bPCA <- function(loadings.chain, 
                                  vars.to.get, 
                                  axes.to.get,
                                  quantiles=c(0.025, 0.5, 0.975))
{
  require(reshape)
  require(lattice)
  V <- nrow(loadings.chain[[1]])
  
  # if vars.to.plot not specified, plotting all (but max 5) variables
  if(missing(vars.to.get)) vars.to.get <- 1:min(V, 5)
  if(missing(axes.to.get)) axes.to.get <- 1:min(V, 5)
  
  # function trimmer() that deletes undesired vars and axes
  trimmer <- function(M, vars.to.get, axes.to.get){
    M[vars.to.get, axes.to.get]
  }
  # apply the trimmer() to the chain of loadings (loop avoidance)
  loadings.trimmed <- lapply(loadings.chain, 
                             FUN=trimmer, vars.to.get, axes.to.get) 
  
  # "melt" the chain into a data frame
  melted <- melt(loadings.trimmed)[,1:3]
  
  # plot using lattice graphics
  print(histogram(~ value | X2 * X1, data=melted,
                  panel = function(x, ...){
                    panel.histogram(x, ...)
                  }))
  
  # summarizing the data for output using quantiles
  output <- vector(mode="list", length=length(quantiles))
  names(output) <- as.character(quantiles)
  for(i in 1:length(quantiles))
  {
    output[[i]] <- apply(simplify2array(loadings.trimmed), 1:2, 
                         FUN=quantile, 
                         probs=quantiles[i])
  }
  return(output)
}

# ------------------------------------------------------------------------------
# FAST LOOP-AVOIDING FUNCTION THAT EXTRACTS THE CHAINS FOR THE SCORES

get.scores.chain.bPCA <- function(bPCA.fitted, data)
{
  sims <- bPCA.fitted$BUGSoutput$sims.matrix
  V = ncol(data)
  
  # split by rows (for the lapply function below)
  sims <- split(sims, seq(nrow(sims)))
  names(sims) <- NULL
  
  # the function that will be used by lapply
  scores.extract <- function(sim, V, data)
  {
    covm = matrix(sim[1:(V^2)], V, V)
    mus <- sim[((V^2)+1):((V^2)+V)]
    centered = scale(data, center=mus, scale=FALSE)
    loadings = eigen(covm)$vectors
    scores = centered  %*% loadings
    
    row.names(scores) <- row.names(data)
    colnames(scores) <- paste("Comp.", 1:V, sep="")
    return(scores)
  }
  
  scores.chain <- lapply(X=sims, FUN=scores.extract, V=V, data=data)
  return(scores.chain)  
}

# ------------------------------------------------------------------------------
# THIS FUNCTION SUMMARIZES THE CHAIN OF THE SCORES BY QUANTILES, 
# NOTHING IS PLOTTED.

summary.scores.bPCA <- function(scores.chain, axes.to.get, 
                                quantiles=c(0.025, 0.5, 0.975))
{
  # function trimmer2() that deletes undesired columns (axes) in a matrix:
  trimmer2 <- function(M, axes.to.get){
    M2 <- as.matrix(M[, axes.to.get])
    colnames(M2) <- paste("Comp.", axes.to.get, sep="")
    return(M2)
  }
  # apply the trimmer() to the chain of loadings (loop avoidance)
  scores.trimmed <- lapply(scores.chain, 
                           FUN=trimmer2, axes.to.get)
  
  #summarizing the data for output using quantiles
  output <- vector(mode="list", length=length(quantiles))
  names(output) <- as.character(quantiles)
  for(i in 1:length(quantiles))
  {
    # calculating the quantile - again, loop avoidance
    output[[i]] <- apply(simplify2array(scores.trimmed), 1:2, 
                         FUN=quantile, 
                         probs=quantiles[i])
  }
  return(output)
}

# summarizing the first axis
