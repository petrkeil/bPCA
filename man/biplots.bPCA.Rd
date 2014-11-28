\name{biplots.bPCA}
\alias{biplots.bPCA}
\title{
bPCA biplots
}

\arguments{
  \item{bPCA.fitted}{
Object of class \code{rjags} returned by the \code{sim.bPCA} function.
}
  \item{data}{
Data matrix in standard form for ordination analyses.
}
  \item{axes.to.plot}{
Vector of indices of principal components that are plotted.
}
  \item{scale}{
Scaling constant.
}

}

\details{
The function returns three PCA biplots. Two are for the lower and upper 95 percentiles, and one is for median of the MCMC simulations. The idea is to show how stable is the PCA.
}

\usage{
biplots.bPCA(bPCA.fitted, data, axes.to.plot=1:2, scale=1)
}

\examples{
bPCA.fitted <- sim.bPCA(iris[,1:4])
biplots.bPCA(bPCA.fitted, iris[,1:4])
}