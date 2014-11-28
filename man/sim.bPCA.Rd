\name{sim.bPCA}
\alias{sim.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function that fits multivariate normal distribution to the data
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
Data matrix in standard form for ordination analyses.
}
  \item{covmat.prior}{
Covariance matrix prior (note that whole covariance matrix must be provided, not single covariances).
}
  \item{covmat.prior.DF}{
Informativeness of covariance matrix prior (integer).
}
  \item{mu.prior}{
Vector of prior means of multivariate normal distribution.
}
  \item{mu.prior.cov}{
mu prior informativeness (matrix of covariances).
}
  \item{n.chains}{
Number of chains in the MCMC algorithm.
}
  \item{n.iter}{
Total number of iterations in the MCMC algorithm.
}
  \item{n.burnin}{
Number of burn-in iterations in the MCMC algorithm.
}
}
\details{
This is the most important function in the \code{bPCA} package. The \code{sim.bPCA} function runs MCMC sampling to fit multivariate normal distribution to the data. The function uses a simple JAGS script, and may take some time to run over large datasets.

This is the BUGS code of the model:

\code{
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
}
  
}
\value{
\code{sim.bPCA} returns an \code{rjags} object which contains MCMC chains for the vector of MVN (Multi-Variate Normal) means and MVN covariance matrix.
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
bPCA.fitted <- sim.bPCA(iris[,1:4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
