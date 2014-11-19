\name{summary.scores.bPCA}
\alias{summary.scores.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes the postriors of the bPCA scores by quantiles
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.scores.bPCA(scores.chain, axes.to.get, quantiles = c(0.025, 0.5, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores.chain}{
%%     ~~Describe \code{scores.chain} here~~
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
function (scores.chain, axes.to.get, quantiles = c(0.025, 0.5, 
    0.975)) 
{
    trimmer2 <- function(M, axes.to.get) {
        M2 <- as.matrix(M[, axes.to.get])
        colnames(M2) <- paste("Comp.", axes.to.get, sep = "")
        return(M2)
    }
    scores.trimmed <- lapply(scores.chain, FUN = trimmer2, axes.to.get)
    output <- vector(mode = "list", length = length(quantiles))
    names(output) <- as.character(quantiles)
    for (i in 1:length(quantiles)) {
        output[[i]] <- apply(simplify2array(scores.trimmed), 
            1:2, FUN = quantile, probs = quantiles[i])
    }
    return(output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
