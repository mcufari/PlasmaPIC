#ifndef GRID_HEADER
#define GRID_HEADER
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "../particles/Particle.hpp"
#include "mkl.h"
#include <immintrin.h>

template<int dim = 0>
class Grid{
public:
    Grid(); 
};


template<> class Grid<1>{
    public:
    Grid(int NPoints, double dx, double lBoundary);
 
    //Grid(int NPoints, double dx, double lBoundary, std::string boundCondition);
  
    void vertexInfoTraverse();
    void poissonSolver();
    void EYBZSolver();

    void getInitialVelocities(double* pList, int* iList, int NParticles);
    void updateVelocities(double* pList, int* iList, int NParticles);
    void moveParticles(double* pList, int* iList, int NParticles);
    void particleInCell(double* pList, int* indexArr, const int NParts);
    void currentInCell(double* pList, int* iList, const int NParts, double* currentList);

    void particleInitRandomStaticProton(Particle<1>* pList, const int NParts) const;
    void IntegrationLoop(double* pList, int* iList, const int NParts);
    void Initialize(double* pList, const int NParts, int* iList);
    void particleInfoTraverse(double* pList, int* iList, int NParticles);
    void particleInitStaticProtonPair(Particle<1>* pList, const int NParts);
    void particleInitUniformProtonElectronPairs(Particle<1>* pList, const int NParts);
    void particleInitSinusoidElectrons(double* pList, int* iList, const int NParts);

    double* jyUpValues;
    double* jyDownValues;
    // double* jzUpValues;
    // double* jzDownValues;

    double* fzUpValues;
    double* fzDownValues;
    double* gridLocations;
    double* EfieldValues;
    double* EFieldY;
    double* EFieldZ;

    double* BfieldValues;
    double* BFieldY;
    double* BFieldZ;
    double* chargeDensity;    

    int NP;
    double dx;
    double lBound;
    double rBound;
    double dt;
    double* poissonDiagonal;
    double* poissonUpDiagonal;
    double* poissonLDiagonal;
    double* phi;
    int* ipiv;
    std::string boundaryCondition;
    long int timestep;
    int timeStepRate = 100000;
    int dump = 1;
    double time = 0.0;
};


#endif 