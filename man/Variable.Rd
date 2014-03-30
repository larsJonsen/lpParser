\name{Variable}
\alias{Variable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Define variables of the model
}
\description{
Variables are the unknown values that are fund by solving the linear programming model. They have dimension according to the indexes
}
\usage{
Variable(model) <- list
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A lpParser object with defined index.
  }
  \item{list}{
  A named list with \code{varable-name = c(index)} with index being a vector of index.
  }
}
\details{
\code{Variable} modifies the lpParser object. Each time \code{Variable()} is called it will add more variables to the model. There is know way to remove them.}
\value{
\code{Variable} is a replacement methods that modify the original object. The function will add to the variable slots in the model. No values are returned.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lars Ravn-Jonsen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Index}}, \code{\link{Variable}}, \code{\link{Constant}}, \code{\link{Equation}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
model <- lpParser()

Index(model) <- list(time = 1:1000, tecnology = c('A', 'B', 'C', 'D'))

Variable(model) <- list(production = c('time', 'tecnology'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{variablde}
\keyword{linear programing}% __ONLY ONE__ keyword per line
