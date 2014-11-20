\name{get.loadings.chain.bPCA}
\alias{get.loadings.chain.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts chains of the loagings from the fitted bPCA object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.loadings.chain.bPCA(bPCA.fitted, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bPCA.fitted}{
object of class \code{rjags} returned by the \code{sim.bPCA} function.
}
  \item{data}{
data matrix in standard form for ordination analyses
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
function (bPCA.fitted, data) 
{
    sims <- bPCA.fitted$BUGSoutput$sims.matrix
    V <- ncol(data)
    sims.cov <- sims[, 1:(V^2)]
    sims.cov <- split(sims.cov, seq(nrow(sims.cov)))
    names(sims.cov) <- NULL
    load.extract <- function(cov, V, data) {
        covm = matrix(cov, V, V)
        loadings <- eigen(covm)$vectors
        row.names(loadings) <- names(data)
        colnames(loadings) <- paste("Comp.", 1:V, sep = "")
        return(loadings)
    }
    loadings.chain <- lapply(X = sims.cov, FUN = load.extract, 
        V = V, data = data)
    return(loadings.chain)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
