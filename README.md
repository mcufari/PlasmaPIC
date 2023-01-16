# PlasmaPIC

A work in progress particle-in-cell plasma simulation code

PlasmaPIC is (planned) multidimensional, capable of simulating a wide array of problems in various geometries

## PIC 1-D

PIC 1-D is the 1d3v simulation module of plasmaPIC. The 1D poisson-solver uses LAPACK to solve the system of equations from the second-order finite difference scheme, then performs finite differencing of the potential to obtain the electric field. The EY and BZ field components are calculated from the Jy current density

## PIC 2-D (Not yet implemented)

PIC 2-D is the 2-dimensional simulation module of plasma PIC. There are two primary options for poisson solvers: FFT (provided from intel MKL) and intel MKL poisson solver. 

## PIC 3-D (Not yet implemented)

PIC 3-D is the 3-dimensional simulation module of plasma PIC. There are again two primary options for poisson solvers: FFT and intel MKL poisson solver. Either can be specified in the input/setup routine


## Test-suite (In progress)

Verify the performance and accuracy of PlasmaPIC on a number of test-problems
  1. Two +1 charge particles
  2. Sole particle (No self-forces condition)
  3. Uniform cold plasma (electron oscillations)


### Compiling and running plasmaPIC

PlasmaPIC has been tested on MacOS with the latest iteration of intel MKL using the intel classic C compiler (icc).
Clone the repository, source oneapi/setvars.sh, use the command  "make plasma" in the repository.

PlasmaPIC will use cmake in the future



### Planned features
  Multithreading with OpenMP and MPI for cluster computing environments
  GPU offloading for speedups of poisson/interpolation bottlenecks

PlasmaPIC has not yet been tested with OpenMP^*
PlasmaPIC has not yet been tested with GPU offloading


