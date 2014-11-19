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
%%     ~~Describe \code{loadings.chain} here~~
}
  \item{vars.to.get}{
%%     ~~Describe \code{vars.to.get} here~~
}
  \item{axes.to.get}{
%%     ~~Describe \code{axes.to.get} here~~
}
  \item{quantiles}{
%%     ~~Describe \code{quantiles} here~~
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
function (loadings.chain, vars.to.get, axes.to.get, quantiles = c(0.025, 
    0.5, 0.975)) 
{
    require(reshape)
    require(lattice)
    V <- nrow(loadings.chain[[1]])
    if (missing(vars.to.get)) 
        vars.to.get <- 1:min(V, 5)
    if (missing(axes.to.get)) 
        axes.to.get <- 1:min(V, 5)
    trimmer <- function(M, vars.to.get, axes.to.get) {
        M[vars.to.get, axes.to.get]
    }
    loadings.trimmed <- lapply(loadings.chain, FUN = trimmer, 
        vars.to.get, axes.to.get)
    melted <- melt(loadings.trimmed)[, 1:3]
    print(histogram(~value | X2 * X1, data = melted, panel = function(x, 
        ...) {
        panel.histogram(x, ...)
    }))
    output <- vector(mode = "list", length = length(quantiles))
    names(output) <- as.character(quantiles)
    for (i in 1:length(quantiles)) {
        output[[i]] <- apply(simplify2array(loadings.trimmed), 
            1:2, FUN = quantile, probs = quantiles[i])
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
