\name{lpParser-package}
\alias{lpParser-package}
\docType{package}
\title{
Interface to linear programming models.
}
\description{
Interface to linear programming models. The package supply functions to parse a optimizing problem written as text string into a linear programming model and solve it with lpSolve via lpSolveAPI.
}
\details{
\tabular{ll}{
Package: \tab lpParser\cr
Type: \tab Package\cr
Version: \tab 0.1.1\cr
Date: \tab 2014-02-08\cr
License: \tab GPL\cr
}

The package gives an easy interface to the lpSolveAPI package for construction of linear programming models based on indexed variable and constants:

Index: \eqn{time : 1:T} and \eqn{technology : A, B, C, D}

Variables: \eqn{production(time, technology)}

Objective: \eqn{\min_{production} (cost production)}{min(cost production}

Subjected to:

\eqn{Demand(time): sum(production) >=  demand}

\eqn{Capacity(time, technology): production <=  capacity}

 \eqn{RampUp(time(1:(T-1)), technology):  production(t - 1,...) - production(t, ...)) <= rampUp(...)}
 
 \eqn{RampDown(time(1:(T-1)): technology)  production(t - 1,...) - production(t, ...)) >= rampDown(...)} 




To simplify the problem notice that a important feature of the model is the dimension of the index; here time and technology. Accordingly this the first to  define. Most of the index are not needed in the constraints (the ...) because they are given by the dimension of equation. Eg the Demand(time) constraint: Here there will be one line for each element of the time dimension, production that have the dimension (time, technology) then have to have one element per technology per constraint.
}
\author{
Lars Ravn-Jonsen

Maintainer: Lars Ravn-Jonsen <lrj@sam.sdu.dk>
}
\references{
~~ Literature or other references for background information ~~
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation ~~
%%~~ directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
