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
Object of class \code{rjags} returned by the \code{sim.bPCA} function.
}
  \item{data}{
Data matrix in standard form for ordination analyses.
}
}
\details{
The returned object can be further summarized by the \code{summary.scores.bPCA} function.
}
\value{
Returns a list of length that is equal to the number of MCMC steps. Each element
is a matrix of the PCA scores.
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
get.scores.chain.bPCA(bPCA.fitted, iris[,1:4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
