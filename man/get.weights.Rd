\name{get.weights}
\alias{get.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract propensity score weights
}
\description{
Extracts propensity score weights from a \code{ps.cont} object.
}
\usage{
get.weights(ps1, stop.method = "wcor", withSampW = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps1}{
a \code{ps.cont} object
}
  \item{stop.method}{
indicates which set of weights to retrieve from the \code{ps.cont} object
}
  \item{withSampW}{
Returns weights with sample weights multiplied in, if they were provided in the original \code{ps.cont} call.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
a vector of weights
}
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
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
