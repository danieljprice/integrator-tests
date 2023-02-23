integrator-tests: a simple code for testing integrators
=======================================================
This code is a simple two-body integrator used to try out various time integration algorithms

Requirements
------------
- Gnu Fortran compiler (part of gcc), e.g. "brew install gcc" or "sudo apt-get install gcc"

Usage:
------
```
make
./nbody
```
the outputs are written in various files called *.out (x and y positions as function of time) and *.ev (energy vs time)
where the name of each file corresponds to the integration method being used. 

The integrators tested are:

- Second order Runge Kutta method
- Second order symplectic Leapfrog method
- Second order leapfrog in velocity verlet form (symplectic)
- Second order verlet with adaptive timesteps
- Fourth order Runge Kutta method
- Forest-Ruth 4th order symplectic
- Position-Extended Forest-Ruth Like (PEFRL) 4th order symplectic

References:
-----------
Young (2012) "The leapfrog method and other ``symplectic'' algorithms for integrating Newton's laws of motion", online lecture notes:

https://young.physics.ucsc.edu/115/leapfrog.pdf

Credits:
--------
Written by Daniel Price <daniel.price@monash.edu>
Code is released under MIT License
