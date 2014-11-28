\name{summary.loadings.bPCA}
\alias{summary.loadings.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes the posteriors of the bPCA loaginds by quantiles and histograms
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.loadings.bPCA(loadings.chain, vars.to.get, axes.to.get, quantiles = c(0.025, 0.5, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loadings.chain}{
Output of function \code{get.loadings.chain.bPCA}
}
  \item{vars.to.get}{
Numeric vector of indices of variables to be summarized.
}
  \item{axes.to.get}{
Numeric vector of indices of axes to be summarizes.
}
  \item{quantiles}{
Numeric vector of quantiles to be used as summaries. The default is 95 percentiles and a median.
}
}
\details{
The function summarizes the MCMC chains of the loadings in two ways:

- It provides a \code{lattice} histograms of all posterior densities of the selected loadings.

- It provides list of quantiles of the selected loadings. 
}
\value{
List of quantiles of the selected loadings. 
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
loadings <- get.loadings.chain.bPCA(bPCA.fitted, iris[,1:4])
summary.loadings.bPCA(loadings)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
