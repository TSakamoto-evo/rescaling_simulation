Code accompanying Sakamoto (2026).
Each folder contains the code used to run simulations with our simulator and SLiM.

Our simulator is written in C++.
In the current version, OpenMP parallelization is used in one part of the simulation code (recombination-mutation event).
To use this version, an OpenMP-compatible compiler is required.
The executable program can be compiled, for example, with:

`g++ *.cpp -Wall -Wextra -std=c++17 -O3 -fopenmp -o XXX.out`

For single-core use, a non-parallel version is also provided and can be compiled with:

`g++ *.cpp -Wall -Wextra -std=c++17 -O3 -o XXX.out`

The only difference between the two versions is the parallelization of a `for` loop in `population.cpp`. 
Therefore, the two versions are intended to have essentially the same simulation behavior.