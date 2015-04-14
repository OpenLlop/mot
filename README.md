MOT v1.0
========

Matlab Optimisation Toolbox (MOT) is a library of optimisation functions for Matlab.

The toolbox includes the following algorithms:

* Genetic Algorithm (**GA**)
* Genetic Algorithm with Islands (**GA** + **Islands**)
* Particle Swarm (**PS**)
* Simulated Annealing (**SA**)

They are intended to be used for optimization of functions of real arguments, integer arguments, a mixture of both or actually any kind of argument, since call-back functions are used. This means that any data type can be used to represent the problem.

The functions are coded focusing on an educational usage; algorithm clarity and ease-of-understanding is given high priority at the expense of any potential code/performance optimisations.

Example scripts are also provided. Some examples show how to improve the result of genetic algorithm with a conventional optimization code.

Code files
----------

| Function | Description | Dependencies
|---------:|:------------|:-----------:
| `aga` | Iterates to find minimum of a function using **GA** | 
| `aga_islands` | Iterates to find minimum of a function using **GA** + **Islands** | `aga`
| `aps` | Iterates to find minimum of a function using **PS** | 
| `asa` | Iterates to find minimum of a function using **SA** | 

Example scripts
---------------

| Script | Description | Dependencies
|--------:|:------------|:-----------:
| `example_aga_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **GA** | `aga`
| `example_aga_islands_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **GA** + **Islands** | `aga` `aga_islands`
| `example_aga_knapsack` | Find local minima of a set of [0,1] integer values (knapsack problem) using **GA** | `aga`
| `example_aps_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **PS** | `aps`
| `example_asa_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **SA** | `asa`

