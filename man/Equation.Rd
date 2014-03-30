\name{Equation}
\alias{Equation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Equation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Equation(model) <- string
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
  A lpParser object with defined index, constants and variables.
}
  \item{string}{
  A sting or a vector of strings each formatted as a equation. See details
}
}
\details{
  A equation \code{string} is a constraint represented by a string with the elements: 
  
  \verb{'name(index); lhs <=> rhs'}

where \code{name} is a name chosen to significant the constraint. \code{Index} is the dimension of the constraint, \code{lhs} is the left hand side,  \code{<=>} is the type of constraint and \code{rhs} is the right hand side. \cr

The left hand side is a linear combination of variables and constants. It consists of \code{constant * variable} separated by a sign (+ or -). The \code{constant} part is optional if missing one is assumed. The \code{constant} can be a declared constant or a number. The type of the constraint can be either '<=', '==' or '=>'. The right hand side can be either a number or a declared constant. The \code{index} of the equation determines the number of equations set up. Constant and variables are expanded according to the index of the equation: If the variable have more dimension than the equation the variable will be represented with one each of that index values in each equation, will the index that corespondents with the equation will be represented with the value of the equation.

Optional the index can be sliced by a [...] statement eg:

\verb{'name(index[1:(t-1)]); constant * variable[1:(t-1)] - constant * variable[2:t] <=> rhs'}

allowing for constraint involving variables with eg time lag.
  }

\value{

Equation is a replacement method that modify the original object. The function will add to the equation slots in the model. No values are returned.
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
\code{\link{Index}}, \code{\link{Variable}}, \code{\link{Constant}}, \code{\link{Solve}}, \code{\link{Sense}}, \code{\link{Bounds}}, \code{\link{lhsMatrix}}, \code{\link{lhsImage}}
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
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
