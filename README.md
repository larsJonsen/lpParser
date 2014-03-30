lpParser
==========

Interface to linear programming models in R. The package supply functions to parse a optimizing problem written as text string into a linear programming model and solve it with lpSolve via `lpSolveAPI`. 

The library is created as an easy interface for learning linear programming model for for example time series production models.

The parser works with the following concepts:

- Index: The index n the equations, for example time: `1:T` or technology: `A, B, C, D`.
- Variables: The variables which values are to be found in the model, They will normally have index. For example `production(time, technology)`.
- Constants: The constants used in the formulation. They can have a index for example `cost(technology)` meaning there is one cost corresponding to the technologies. 
- Objective: The objective to be optimized, it will have to be a linear formulation of the variables. For example `cost * production`.
- Equations: Equations are the constraints for optimisation. They have a left hand side, a type and a right hand side. Equations have name and index. For example `Capacity(time, technology): production <=  capacity` meaning that the production af each technology have to be with in capacity at all time. It will expand into `length(time) * length(technology)` equations:


```
model <- lpParser()

n=10 #reduced model in time dimension
Index(model) <- list(time = 1:n, technology = c('A', 'B', 'C', 'D'))

Constant(model) <- list(
      demand = c('time'), 
      capacity = c('technology'),
      cost = 'technology',
      rampUp = 'technology',
      rampDown ='technology'
)

setConstant(model) <- list(demand = demand[1:n], 
     capacity = powerPlants$Capacity,
     cost = rowSums(powerPlants[,c('opmain', 'fuelcost')]),
     rampUp = powerPlants$Rampup,
     rampDown = powerPlants$Rampdown)
     
Variable(model) <- list(production = c('time', 'technology'))

Equation(model) <- c('Demand(time); production >= demand')
Equation(model) <- c('Capacity(time,technology);  production <= capacity')
Equation(model) <- c('RampUp(time[1:(n-1)],technology); production[2:n] - production[1:(n-1)] <= rampUp')
Equation(model) <- c('RampDown(time[1:(n-1)],technology); production[2:n] - production[1:(n-1)] >= rampDown')

Object(model) <- 'cost * production'

Solve(model)
```
