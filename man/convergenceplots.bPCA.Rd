\name{convergenceplots.bPCA}
\alias{convergenceplots.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
bPCA convergence plots
}
\description{
Simple convergence diagnostics plots calculated from the fitted \code{rjags} object returned by the \code{sim.bPCA()} function.
}
\usage{
convergenceplots.bPCA(bPCA.fitted)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bPCA.fitted}{
object of class \code{rjags} returned by the \code{sim.bPCA} function.
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
bPCA.fitted <- sim.bPCA(iris[,1:4])
convergenceplots.bPCA(bPCA.fitted)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
