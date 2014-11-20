\name{eigenvalplots.bPCA}
\alias{eigenvalplots.bPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that summarizes the bPCA posteriors of eigenvalues and explained variance
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
    V <- ncol(data)
    sims <- bPCA.fitted$BUGSoutput$sims.matrix
    sims <- sims[, 1:(V * V)]
    eigen.chains <- matrix(nrow = nrow(sims), ncol = V)
    for (i in 1:nrow(sims)) {
        covm <- matrix(sims[i, ], V, V)
        eigen.chains[i, ] <- eigen(covm)$values
    }
    exp.vars <- eigen.chains/rowSums(eigen.chains) * 100
    par(mfrow = c(1, 2))
    boxplot(eigen.chains, ylab = "Eigenvalue", xlab = "PCA axis", 
        col = "grey", outline = FALSE)
    boxplot(exp.vars, ylab = "Explained variability [\% of total]", 
        xlab = "PCA axis", col = "grey", outline = FALSE, ylim = c(0, 
            100))
    results <- list(Eigenvalues = summary(eigen.chains), Exp.var = summary(exp.vars))
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
