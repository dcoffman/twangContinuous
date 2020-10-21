\name{summary.ps.cont}
\alias{summary.ps.cont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a ps.cont object
}
\description{
Computes summary information about a stored \code{ps.cont} object
}
\usage{
\method{summary}{ps.cont}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
a \code{\link{ps.cont}} object
}
  \item{\dots}{
additional arguments affecting the summary produced
}
}
\details{
Compresses the information in the \code{desc} component of the \code{ps.cont} object
into a short summary table describing the size of the dataset and the quality of
the propensity score weights.
}
\value{
See \code{\link{ps.cont}} for details on the returned table
}

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
