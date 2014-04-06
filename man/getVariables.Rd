\name{getVariables}
\alias{getVariables}
\alias{getConstraints}
\alias{getObjective}
\title{
Get the solved model
}
\description{
getVariables get the variables form a solved model, getConstraints get the value and shadow values of the constraints from a solved model and getObjective returns the objective from a solved model.
}
\usage{
getVariables(model)
getConstraints(model)
getObjective(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A solved lpParser model
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{getVariable} returns a named list of arrays. The names correspond to the variables names and the dimension of the arrays corresponds to the dissing of the respective variables.

\code{getVariable} returns a name list of data.frames. The names corresponds to the constraints  names and the data frame have the following properties: row names corresponds to the number of the constraint, one column for each of the indices of the constraint named according to the constraint, one column named value giving the value of the constraint, and one column called lambda giving the shadow value of the constraint (the value of the constraint in the dual solution)

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
\code{\link{Index}}, \code{\link{Variable}}, \code{\link{Constant}}, \code{\link{Equation}}, \code{\link{Sense}}, \code{\link{Bounds}}, \code{\link{Solve}}}
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

getVariables(model)
getConstraints(model)
getObjective(model)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
