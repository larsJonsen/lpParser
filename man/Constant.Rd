\name{Constant}
\alias{Constant}
\alias{setConstant}
\title{
Define and set constants
}
\description{
Constants are the known parameters of the model. They have dimensions according to indexes. Constants are first defined according to index and then set with values. If a constant is single valued it have the index \code{'single'}.
}
\usage{
Constant(model) <- nameList
setConstant(model) <- valueList
}
\arguments{
  \item{model}{
  A lpParser object with defined index.
  }
  \item{nameList}{
  A named list with constant-name = c(index) with index being a vector of index.
  }
  \item{valueList}{
  A named list with the constant-names = array with array having the values and with     dimension corresponding to the length of the index.
  }
}

\details{
Both \code{Constant} and \code{setConstans} modifyes the lpParser object. Each time constant is called it will add more constants to the model. There is know way to remove them.
}
\value{
\code{Constant} and \code{setConstant} are replacement methods that modify the original object. The function will add to the constant slots in the model. No values are returned.
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
Lars Ravn-Jonsen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{\code{\link{lpParser-package}} \code{\link{Index}}
}

\examples{
model <- lpParser()

n = 10 # reduced time dimension
Index(model) <- list(time = 1:n, tecnology = c('A', 'B', 'C', 'D'))

Constant(model) <- list(
      demand = c('time'), 
      capacity = c('tecnology'),
      cost = 'tecnology',
      rampUp = 'tecnology',
      rampDown ='tecnology'
)

setConstant(model) <- list(demand = demand[1:n], 
     capacity = powerPlants$Capacity,
     cost = rowSums(powerPlants[,c('opmain', 'fuelcost')]),
     rampUp = powerPlants$Rampup,
     rampDown = powerPlants$Rampdown)
  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{constant}
\keyword{parameters}
\keyword{linear programing}
