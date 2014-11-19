\name{get.scores.chain.bPCA}
\alias{get.scores.chain.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extracts chain of PCA scores from fitted bPCA object
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get.scores.chain.bPCA(bPCA.fitted, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bPCA.fitted}{
%%     ~~Describe \code{bPCA.fitted} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
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
    V = ncol(data)
    sims <- split(sims, seq(nrow(sims)))
    names(sims) <- NULL
    scores.extract <- function(sim, V, data) {
        covm = matrix(sim[1:(V^2)], V, V)
        mus <- sim[((V^2) + 1):((V^2) + V)]
        centered = scale(data, center = mus, scale = FALSE)
        loadings = eigen(covm)$vectors
        scores = centered \%*\% loadings
        row.names(scores) <- row.names(data)
        colnames(scores) <- paste("Comp.", 1:V, sep = "")
        return(scores)
    }
    scores.chain <- lapply(X = sims, FUN = scores.extract, V = V, 
        data = data)
    return(scores.chain)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
