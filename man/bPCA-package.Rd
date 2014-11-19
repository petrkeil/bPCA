\name{bPCA-package}
\alias{bPCA-package}
\alias{bPCA}
\docType{package}
\title{
Bayesian Principal Components Analysis (bPCA)
}
\description{
Allows to calculate PCA in Bayesian setting (using MCMC sampler in JAGS). Prior information on the associations (covariances) between variables can be provided. Calculates full posterior distributions of eigenvalues, explained variances, PCA loadings and PCA scores. Some basic plotting and summary functions are also provided.
}
\details{
\tabular{ll}{
Package: \tab bPCA\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-11-03\cr
License: \tab GNU GPL (>=2)\cr
}
STEP 1:

The bPCA analysis always starts with a single dataset (\code{data.frame} or \code{matrix}). Each column must be numerical variable.

STEP 2:

The data are first subjected to the \code{sim.bPCA} function which runs MCMC sampling to fit multivariate normal distribution to the data. The function uses a simple JAGS script, and may take some time to run over large datasets. \code{sim.bPCA} returns an \code{rjags} obects which contains MCMC chains for the vector of MVN (Multi-Variate Normal) means and MVN covariance matrix.

STEP 3: 

Various aspects of the Principal Components Analysis can be extracted from the \code{rjags} object. These are:

- Posterior distributions of eigenvalues and explained variances using the \code{eigenvalplots.bPCA}.

- Posterior bPCA biplots using the \code{biplots.bPCA} function.

- Posterior distributions of the loadings using the \code{get.loadings.chain.bPCA} and \code{summary.loadings.bPCA}.

- Posterior distributions of the PCA scores using the \code{get.scores.chain.bPCA} and \code{summary.scores.bPCA}.

}
\author{
Jan Smycka, Petr Keil

Maintainer: Petr Keil <pkeil@seznam.cz>, Jan Smycka <smyckaj@gmail.com>
}
\references{
Legendre & Legendre (1998) Numerical Ecology. Elsevier, Amsterdam.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R ~~
~~ documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{

# We will use first four variables the "iris" dataset in the package "datasets":
summary(iris[,1:4])
head(iris[,1:4])

# Running the MCMC simulation and saving the chains as "bPCA.fitted":
bPCA.fitted <- sim.bPCA(iris[,1:4])

# Checking convergence of the MCMC chains:
convergenceplots.bPCA(bPCA.fitted)

# Posterior distributions of eigenvalues and percentages of
# explained variance:
eigenvalplots.bPCA(bPCA.fitted, iris[,1:4])

# PCA biplots of the 1st and 2nd axes. We extract 2.5%, 50% and 97.5%
# quantiles of the posterior distributions:
biplots.bPCA(bPCA.fitted, iris[,1:4], axes.to.plot=c(1,2), ylim=c(-2,2))

# We can compare the previous plots with classical PCA biplot:
plot.classicPCA(iris[,1:4], axes.to.plot=1:2, scale=1, ylim=c(-1,1))

# Extracting MCMC chains for PCA loadings from the bPCA.fitted object:
loadings.chain <- get.loadings.chain.bPCA(bPCA.fitted, iris[,1:4])

# Quantiles and summary plots of the posteriors of the loadings:
summary.loadings.bPCA(loadings.chain, vars.to.get=1:4, axes.to.get=1:2)

# Extracting MCMC chains for PCA scores from the bPCA.fitted object:
scores.chain <- get.scores.chain.bPCA(bPCA.fitted, iris[,1:4])

# Quantiles of the posteriors of the scores:
summary.scores.bPCA(scores.chain, axes.to.get=1:2)

}
