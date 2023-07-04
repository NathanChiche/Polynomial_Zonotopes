# Polynomial_Zonotopes

One has to download Julia for the code to be executed and the following packages: Nemo, LazySets, Plots, Random, Symbolics, Dates, BenchmarkTools, TimerOutputs

The code runs as a notebook: you compile every package and function and do you main tests in the main function.

One can also include the file in the Julia Repl (might not work at the moment due to type problems)

The goal is to compute an invariant for a given discrete dynamical system. Hence for executing, you should only modify the main function by adding iterate_polynomials_over_PZ([p1,p2],P1,4,0,R,150,20000,12000,1.1,false) where [p1,p2] is the list of polynomials that describe your system.
P1 is the set describing the initial conditions (to modify as you want), 4 is the number of iterations to find the invariant. Next parameters will be explained later
