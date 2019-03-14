\name{Solve}
\alias{Solve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Solve
}
\description{
Solves a lpParser object
}
\usage{
Solve(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A lpParser object with defined index, constants, variables and equations.
}
}
\details{
When index, constants, variables and equations are set then the model can be solved by lpSolve by calling \code{Solve()}. Note that if Sense is not called the model is assumed to be minimized and if Bounds is not called all variables is assumed \eqn{latex}{\ge 0}
}
\value{
Solved modifies the original object if successfully solved. The function will add to the slots in the model. As a security the model is returned invisibly; a standard replacement of the model will work as well.
The status of the model is printed.
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
\code{\link{Index}}, \code{\link{Variable}}, \code{\link{Constant}}, \code{\link{Equation}}, \code{\link{Sense}}, \code{\link{Bounds}}, \code{\link{getVariables}}, \code{\link{getConstraints}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
model <- lpParser()

n=10 #reduced model in time dimension
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
     
Variable(model) <- list(production = c('time', 'tecnology'))

Equation(model) <- c('Demand(time); production >= demand')
Equation(model) <- c('Capacity(time,tecnology);  production <= capacity')
Equation(model) <- c('RampUp(time[1:(n-1)],tecnology); production[2:n] - production[1:(n-1)] <= rampUp')
Equation(model) <- c('RampDown(time[1:(n-1)],tecnology); production[2:n] - production[1:(n-1)] >= rampDown')

Object(model) <- 'cost * production'

Solve(model)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
