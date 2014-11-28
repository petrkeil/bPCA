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
output of function get.scores.chain.bPCA
}
  \item{axes.to.get}{
numeric vector of axes to summarize
}
  \item{quantiles}{
numeric vector of quantiles to summarize
}
}
\details{
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
scores <- get.scores.chain.bPCA(bPCA.fitted, iris[,1:4])
summary.scores.bPCA(scores, axes.to.get=1:2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
