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
    timestep = 0;
    phi = (double*) calloc(NP, sizeof(double));
    dt = 0.1 * dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
        EfieldValues[i] = 0;  
        BfieldValues[i] = 0;
        chargeDensity[i] = 0;
    }

    poissonMatrix = (double*) calloc((NPoints-2)*(NPoints-2), sizeof(double));
    for(int i = 1; i < NPoints-3; i++){
        poissonMatrix[i*(NPoints-2) + i] = -2;
        poissonMatrix[i*(NPoints-2) + i+1] = 1;
        poissonMatrix[i*(NPoints-2) + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    

    poissonMatrix[((NPoints-2)*(NPoints-2))-1] = -2;
    poissonMatrix[((NPoints-2)*(NPoints-2))-2] = 1;
    
    for(int i = 0; i < NPoints - 2; i++){
        for(int j = 0; j < NPoints - 2;j++){
            std::cout << poissonMatrix[i*(NPoints-2) + j] << " ";
        }
        std::cout << std::endl;
    }
    ipiv = (int*) calloc(NPoints-2, sizeof(int));
    int info = 0;
    int nrhs = 1;
    int NPl2 = NP-2;
    dgetrf(&NPl2,&NPl2,poissonMatrix,&NPl2,ipiv,&info);
    

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
    timestep = 0;
    phi = (double*) calloc(NP, sizeof(double));
    dt = 0.1*dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
        EfieldValues[i] = 0;  
        BfieldValues[i] = 0;
        chargeDensity[i] = 0;
    }

    poissonMatrix = (double*) calloc((NPoints-2)*(NPoints-2), sizeof(double));
    for(int i = 1; i < NPoints-3; i++){
        poissonMatrix[i*(NPoints-2) + i] = -2;
        poissonMatrix[i*(NPoints-2) + i+1] = 1;
        poissonMatrix[i*(NPoints-2) + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    

    poissonMatrix[((NPoints-2)*(NPoints-2))-1] = -2;
    poissonMatrix[((NPoints-2)*(NPoints-2))-2] = 1;
    
    for(int i = 0; i < NPoints - 2; i++){
        for(int j = 0; j < NPoints - 2;j++){
            std::cout << poissonMatrix[i*(NPoints-2) + j] << " ";
        }
        std::cout << std::endl;
    }
    ipiv = (int*) calloc(NPoints-2, sizeof(int));
    int info = 0;
    int nrhs = 1;
    int NPl2 = NP-2;
    dgetrf(&NPl2,&NPl2,poissonMatrix,&NPl2,ipiv,&info);

}

void Grid<1>::vertexInfoTraverse(){
    std::ofstream ofile;
    std::string filename = "dumps/gridQuants" + std::to_string(timestep) + ".txt";
    ofile.open(filename);
    ofile << "gridPosits " << "chargeDensity " << "EFieldValues "  << "phi" << std::endl;
    for(int i = 0; i < NP; i++){
        ofile << gridLocations[i] << " ";
        ofile << chargeDensity[i] << " ";
        ofile << EfieldValues[i] << " ";
        ofile << phi[i] << std::endl;
    }
    ofile.close();
}

void Grid<1>::particleInfoTraverse(Particle<1>* pList, int NParticles){
    std::ofstream ofile;
    std::string filename = "dumps/particleProps" + std::to_string(timestep) + ".txt";
    ofile.open(filename);
    ofile << "posits " << "vels " << std::endl;
    for(int i = 0; i < NParticles; i++){
        ofile << pList[i].position << " ";
        ofile << pList[i].velocity << std::endl;
    }
    ofile.close();
}

void Grid<1>::poissonSolver(){
    const char trans = 'T';
    double* density = (double*) malloc(sizeof(double) * (NP-2));
    for(int i = 0; i < NP-2; i++){
        density[i] = -1*chargeDensity[i+1] * dx * dx;
        //std::cout << "Density vector at " << i << " is: " << density[i] << std::endl;
    }
   
    int info = 0;
    int nrhs = 1;
    int NPl2 = NP-2;

    
    dgetrs(&trans,&NPl2,&nrhs,poissonMatrix,&NPl2,ipiv,density,&NPl2,&info);
    
    
    phi[0] = 0;
    phi[NP-1] = 0;
    for(int i = 0; i < NP-2; i++){
        phi[i+1] = density[i];
    }

    EfieldValues[0] = (phi[NP-1] - phi[1])/(2.0*dx);
    EfieldValues[NP-1] = (phi[NP-2] - phi[0])/(2.0*dx);
    for(int i = 1; i < NP-1; i++){
        EfieldValues[i] = (phi[i-1] - phi[i+1])/(2.0*dx);
    }
    
   
    free(density);
}

void Grid<1>::getInitialVelocities(Particle<1>* pList, int NParticles){
    for(int i = 0; i < NParticles; i++){
        int lIndex = floor((pList[i].position - lBound)/dx);
        int rIndex = (lIndex + 1) % NP;
        double effectiveEfield = EfieldValues[lIndex] * (gridLocations[rIndex] - pList[i].position)/dx + 
                                 EfieldValues[rIndex] * (pList[i].position - gridLocations[lIndex])/dx;
        pList[i].velocity -= (pList[i].charge / pList[i].mass) * effectiveEfield * (dt/2.0);
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
        if(pList[i].position < lBound)
            pList[i].position = rBound + (lBound - pList[i].position);
    }
    timestep++;
}

void Grid<1>::Initialize(Particle<1>* pList, const int NParticles){
    particleInitStaticProtonPair(pList, NParticles);
    
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

        particleInfoTraverse(pList, NParticles);
}