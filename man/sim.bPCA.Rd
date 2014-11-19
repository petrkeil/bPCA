\name{sim.bPCA}
\alias{sim.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function that fits multivariate normal to the data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.bPCA(data, covmat.prior, covmat.prior.DF, mu.prior, mu.prior.cov, n.chains, n.iter, n.burnin, simulation)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{covmat.prior}{
%%     ~~Describe \code{covmat.prior} here~~
}
  \item{covmat.prior.DF}{
%%     ~~Describe \code{covmat.prior.DF} here~~
}
  \item{mu.prior}{
%%     ~~Describe \code{mu.prior} here~~
}
  \item{mu.prior.cov}{
%%     ~~Describe \code{mu.prior.cov} here~~
}
  \item{n.chains}{
%%     ~~Describe \code{n.chains} here~~
}
  \item{n.iter}{
%%     ~~Describe \code{n.iter} here~~
}
  \item{n.burnin}{
%%     ~~Describe \code{n.burnin} here~~
}
  \item{simulation}{
%%     ~~Describe \code{simulation} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, covmat.prior, covmat.prior.DF, mu.prior, mu.prior.cov, 
    n.chains, n.iter, n.burnin, simulation) 
{
    require(R2jags)
    require(MASS)
    require(Matrix)
    require(coda)
    N = nrow(data)
    V = ncol(data)
    if (missing(covmat.prior)) 
        covmat.prior = as.matrix(Diagonal(V, 1/1000))
    if (missing(covmat.prior.DF)) 
        covmat.prior.DF = V
    if (missing(mu.prior)) 
        mu.prior = rep(0, V)
    if (missing(mu.prior.cov)) 
        mu.prior.cov = as.matrix(Diagonal(V, 1000))
    if (missing(n.chains)) 
        n.chains = 3
    if (missing(n.iter)) 
        n.iter = 5000
    if (missing(n.burnin)) 
        n.burnin = 4500
    mu.prior.prec = ginv(mu.prior.cov)
    listdata = list(Y = as.matrix(data), N = N, V = V, covmat.prior = covmat.prior, 
        mu.prior = mu.prior, covmat.prior.DF = covmat.prior.DF, 
        mu.prior.prec = mu.prior.prec)
    cat("\n  model\n  {\n    # priors on the vector of multinormal means\n    mu[1:V] ~ dmnorm(mu.prior[], mu.prior.prec[,])\n    # priors on the covariance matrix\n    prec[1:V,1:V] ~ dwish(covmat.prior[,], covmat.prior.DF)\n    \n    # makes covariance from precision\n    cov[1:V,1:V] <- inverse(prec[,])\n\n    # likelihood\n    for (i in 1:N)\n    {\n      Y[i,1:V] ~ dmnorm(mu[], prec[,])\n    }\n  }\n  ", 
        file = "PCA.bugs")
    pcabay <- jags(data = listdata, model.file = "PCA.bugs", 
        parameters.to.save = c("cov", "mu"), n.chains = n.chains, 
        n.iter = n.iter, n.burnin = n.burnin, DIC = FALSE)
    return(pcabay)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
