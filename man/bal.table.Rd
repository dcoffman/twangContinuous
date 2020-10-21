\name{bal.table}
\alias{bal.table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute balance table
}
\description{
Extract the balance table from \code{\link{ps.cont}} object
}
\usage{
bal.table(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{\link{ps.cont}} object}
  \item{digits}{The number of digits that the numerical entries should be rounded to.}
  \item{\dots}{Additional arguments that may be passed to the underlying function}}
\details{
\code{bal.table} is a generic function for extracting balance
tables from a \code{\link{ps.cont}} object. This object
usually has several sets of candidate weights, one for an unweighted
analysis and perhaps several \code{stop.methods}. \code{bal.table}
will return a table for each set of weights combined into a list. Each list
component will be named as given in the \code{x}, usually the name of the
\code{stop.method}. The balance table labeled \dQuote{unw} indicates the
unweighted analysis.
}
\value{
Returns a data frame containing the balance information.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Donna L. Coffman
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
The example for \code{\link{ps.cont}} contains an example of the
use of \code{bal.table}
}
\examples{
##---- Should be DIRECTLY executable !! ----
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }% use one of  RShowDoc("KEYWORDS")
