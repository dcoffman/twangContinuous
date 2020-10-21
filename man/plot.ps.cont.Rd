\name{plot.ps.cont}
\alias{plot.ps.cont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots for \code{ps.cont} objects
}
\description{
This function produces a collection of diagnostic plots for \code{ps.cont} objects.
}
\usage{
\method{plot}{ps.cont}(x, plots = "optimize", subset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \code{ps.cont} object}
  \item{plots}{An indicator of which type of plot is desired. The options are
    \itemize{
      \item \code{"optimize"} A plot of the balance criteria as a function of the GBM
            iteration. The default.
      \item \code{"es"} Plots of the standardized effect size of the pre-treatment
            variables before and after reweighing}}
  \item{subset}{Used to restrict which of the \code{stop.method}s will be used
                in the figure. For example \code{subset = c(1,3)} would indicate that the
                first and third \code{stop.method}s (in alphabetical order of those specified
                in the original call to the ps.cont function) should be included in the
                figure.}
  \item{\dots}{Additional arguments that may be passed to the underlying \code{lattice}
                package plotting functions}
}
\details{
This function produces lattice-style graphics of diagnostic plots.
}
%%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%%}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
Donna L. Coffman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ps.cont}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
\keyword{ multivariate }% __ONLY ONE__ keyword per line
