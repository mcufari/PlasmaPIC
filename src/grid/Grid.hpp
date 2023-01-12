#ifndef GRID_HEADER
#define GRID_HEADER
#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "../particles/Particle.hpp"
#include "mkl.h"

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

    void getInitialVelocities(Particle<1>* pList, int NParticles);
    void updateVelocities(Particle<1>* pList, int NParticles);
    void moveParticles(Particle<1>* pList, int NParticles);
    void particleInCell(const Particle<1>* pList, const int NParts);
    void currentInCell(const Particle<1>* pList, const int NParts, double* currentList);

    void particleInitRandomStaticProton(Particle<1>* pList, const int NParts) const;
    void IntegrationLoop(Particle<1>* pList, const int NParts);
    void Initialize(Particle<1>* pList, const int NParts);
    void particleInfoTraverse(Particle<1>* pList, int NParticles);
    void particleInitStaticProtonPair(Particle<1>* pList, const int NParts);
    void particleInitUniformProtonElectronPairs(Particle<1>* pList, const int NParts);
    void particleInitSinusoidElectrons(Particle<1>* pList, const int NParts);

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
    int timeStepRate = 50;
    int dump = 1;
    double time = 0.0;
};


#endif 