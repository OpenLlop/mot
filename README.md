MOT v1.0
========

Matlab Optimisation Toolbox (MOT) is a library of optimisation functions for Matlab.

The toolbox includes the following algorithms:

* Genetic algorithm (**GA**)
* Genetic algorithm with islands (**GA** + **Islands**)
* Simulated annealing (**SA**)

They are intended to be used for optimization of functions of real arguments, integer arguments, a mixture of both or actually any kind of argument, since call-back functions are used. This means that any data type can be used to represent the problem.

Example scripts are provided. Some examples show how to improve the result of genetic algorithm with a conventional optimization code.

Code files
----------

| Function | Description
|---------:|:-----------
| `aga` | Iterates to find minimum of a function using **GA**
| `aga_islands` | Iterates to find minimum of a function using **GA** + **Islands**
| `asa` | Iterates to find minimum of a function using **SA**

Example scripts
---------------

| Example | Description
|--------:|:-----------
| `example_aga_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **GA**
| `example_aga_islands_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **GA** + **Islands**
| `example_asa_rastrigin` | Find local minima of a R^2->R function based on _Rastrigin_ function using **SA**

