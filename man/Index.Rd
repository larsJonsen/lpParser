\name{Index}
\alias{Index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replacement method for setting index in the model
}
\description{
Index are the dimensions of the model. It corespondents to the index in the mathematical model. There is a standard index called 'single' that is of length 1.
}
\usage{
Index(model) <- list
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A lpParser object.
  }
  \item{list}{
Index takes as argument a list with index = values where values are the values the index can take.
}

}
\value{
Index is a replacement method that modify the original object. The function will add to the index slots in the model. No values are returned.
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

\seealso{\link{lpParser-package}, \code{\link{Variable}}, \code{\link{Constant}}, \code{\link{Equation}}
}
\examples{
model <- lpParser()

Index(model) <- list(time = 1:1000, tecnology = c('A', 'B', 'C', 'D'))

}
\keyword{index}
\keyword{linear programing}
