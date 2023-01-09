#ifndef GRID_HEADER
#define GRID_HEADER
#include <iostream>
#include <string>
#include "mkl.h"
template<int dim = 0>
class Grid{
public:
    Grid(); 
};


template<> class Grid<1>{
    public:
    Grid(int NPoints, double dx, double lBoundary);
 
    Grid(int NPoints, double dx, double lBoundary, std::string boundCondition);
  
    void vertexInfoTraverse();
    void poissonSolver();
    void getInitialVelocities(Particle<1>* pList);
    
    double* gridLocations;
    double* EfieldValues;
    double* BfieldValues;
    double* chargeDensity;    

    int NP;
    double dx;
    double lBound;
    double rBound;
    double* poissonMatrix;
    int* ipiv;
    std::string boundaryCondition;

    
};


#endif 