/*
Implementations of multidimensional grids
See grid.hpp for documentation
*/

#include "Grid.hpp"
#include <math.h>

Grid<1>::Grid(int NPoints, double dx, double lBoundary) : NP(NPoints), dx(dx), lBound(lBoundary){
    rBound = lBound + (NPoints-1)*dx; //1 point is used at each boundary
    boundaryCondition = "Periodic"; //Only boundary condition permitted at present in periodic
    gridLocations = (double*) calloc(NP, sizeof(double));
    EfieldValues = (double*) calloc(NP, sizeof(double));
    BfieldValues = (double*) calloc(NP, sizeof(double));
    chargeDensity = (double*) calloc(NP, sizeof(double));
    dt = 0.1 * dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
        EfieldValues[i] = 0;  
        BfieldValues[i] = 0;
        chargeDensity[i] = 0;
    }

    poissonMatrix = (double*) calloc(NPoints*NPoints, sizeof(double));
    for(int i = 1; i < NPoints-1; i++){
        poissonMatrix[i*NPoints + i] = -2;
        poissonMatrix[i*NPoints + i+1] = 1;
        poissonMatrix[i*NPoints + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    poissonMatrix[NPoints - 1] = -1;

    poissonMatrix[(NPoints*NPoints)-1] = -2;
    poissonMatrix[(NPoints*NPoints)-2] = 1;
    poissonMatrix[(NPoints-1)*NPoints] = -1;
    
    ipiv = (int*) calloc(NPoints, sizeof(int));
    int info = 0;
    int nrhs = 1;
    
    dgetrf(&NP,&NP,poissonMatrix,&NP,ipiv,&info);
    

}



Grid<1>::Grid(int NPoints, double dx, double lBoundary, std::string boundCondition) : 
NP(NPoints), 
dx(dx), 
lBound(lBoundary), 
boundaryCondition(boundCondition) 
{
    rBound = lBound + (NPoints - 1) * dx;
    gridLocations = (double*) calloc(NP, sizeof(double));
    EfieldValues = (double*) calloc(NP, sizeof(double));
    BfieldValues = (double*) calloc(NP, sizeof(double));
    chargeDensity = (double*) calloc(NP, sizeof(double));
    dt = 0.1*dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
        EfieldValues[i] = 0;  
        BfieldValues[i] = 0;
        chargeDensity[i] = 0;
    }

    poissonMatrix = (double*) calloc(NPoints*NPoints, sizeof(double));
    for(int i = 1; i < NPoints-1; i++){
        poissonMatrix[i*NPoints + i] = -2;
        poissonMatrix[i*NPoints + i+1] = 1;
        poissonMatrix[i*NPoints + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    poissonMatrix[NPoints - 1] = -1;

    poissonMatrix[(NPoints*NPoints)-1] = -2;
    poissonMatrix[(NPoints*NPoints)-2] = 1;
    poissonMatrix[(NPoints-1)*NPoints] = -1;
    
    ipiv = (int*) calloc(NPoints, sizeof(int));
    int info = 0;
    int nrhs = 1;
   
    
     dgetrf(&NP,&NP,poissonMatrix,&NP,ipiv,&info);

}

void Grid<1>::vertexInfoTraverse(){
    std::cout << "gridPosits " << "chargeDensity " << "EFieldValues" << std::endl;
    for(int i = 0; i < NP; i++){
        std::cout << gridLocations[i];
        std::cout << chargeDensity[i];
        std::cout << EfieldValues[i] << std::endl;
    }
}

void Grid<1>::poissonSolver(){
    const char trans = 'T';
    double* density = (double*) malloc(sizeof(double) * NP);
    for(int i = 0; i < NP; i++){
        density[i] = -1*chargeDensity[i] * dx * dx;
        //std::cout << "Density vector at " << i << " is: " << density[i] << std::endl;
    }
   
    int info = 0;
    int nrhs = 1;


    
    dgetrs(&trans,&NP,&nrhs,poissonMatrix,&NP,ipiv,density,&NP,&info);
    
    for(int i = 1; i < NP-1; i++){
        EfieldValues[i] = (density[i-1] - density[i+1])/(2.0*dx);
    }
    EfieldValues[0] = (density[NP-1] - density[1])/(2.0*dx);
    EfieldValues[NP-1] = (density[NP-2] - density[0])/(2.0*dx);
    free(density);
}

void Grid<1>::getInitialVelocities(Particle<1>* pList, int NParticles){
    for(int i = 0; i < NParticles; i++){
        int lIndex = floor((pList[i].position - lBound)/dx);
        int rIndex = (lIndex + 1) % NP;
        double effectiveEfield = EfieldValues[lIndex] * (gridLocations[rIndex] - pList[i].position)/dx + 
                                 EfieldValues[rIndex] * (pList[i].position - gridLocations[lIndex])/dx;
        pList[i].velocity -= (pList[i].charge / pList[i].mass) * effectiveEfield * (dt/2.0);
        std::cout << "Particle " << i << " velocity at -dt/2: " << pList[i].velocity << std::endl; 
    }
}

void Grid<1>::updateVelocities(Particle<1>* pList, int NParticles){
    for(int i = 0; i < NParticles; i++){
        int lIndex = floor((pList[i].position-lBound)/dx);
        int rIndex = (lIndex + 1) % NP;
        double effectiveEfield = EfieldValues[lIndex] * (gridLocations[rIndex] - pList[i].position)/dx + 
                                 EfieldValues[rIndex] * (pList[i].position - gridLocations[lIndex])/dx;
        pList[i].velocity += (pList[i].charge/pList[i].mass) * effectiveEfield * (dt/2.0);
        //Rotate velocities
        pList[i].velocity += (pList[i].charge/pList[i].mass) * effectiveEfield * (dt/2.0);
    }
}

void Grid<1>::moveParticles(Particle<1>* pList, int NParticles){
    for(int i = 0; i < NParticles; i++){
        pList[i].position += pList[i].velocity * dt;
        if(pList[i].position > rBound) 
            pList[i].position = lBound + (pList[i].position-rBound);
        
    }
}

void Grid<1>::Initialize(Particle<1>* pList, const int NParticles){
    particleInitRandomStaticProton(pList, NParticles);
    
    particleInCell(pList, NParticles);

    poissonSolver();

    vertexInfoTraverse();

    getInitialVelocities(pList, NParticles);
}


void Grid<1>::IntegrationLoop(Particle<1>* pList, const int NParticles){
        
        updateVelocities(pList,NParticles);

        moveParticles(pList, NParticles);

        particleInCell(pList, NParticles);

        poissonSolver();

        vertexInfoTraverse();
}