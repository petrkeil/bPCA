\name{eigenvalplots.bPCA}
\alias{eigenvalplots.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that summarizes the bPCA posterior distributions of eigenvalues and explained variances
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
eigenvalplots.bPCA(bPCA.fitted, data)
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
The function gives boxplots that summarize the posterior distributions of eigenvalues and explained variances that correspond with each PC axes.
}
\value{
Apart from the plots the function also returns a list of summary statistics for the eigenvalues and explained variances.
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
  eigenvalplots.bPCA(bPCA.fitted, iris[,1:4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
