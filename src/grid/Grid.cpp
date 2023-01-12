/*
Implementations of multidimensional grids
See grid.hpp for documentation
*/

#include "Grid.hpp"
#include "omp.h"
#include <math.h>

Grid<1>::Grid(int NPoints, double dx, double lBoundary) : NP(NPoints), dx(dx), lBound(lBoundary){
    rBound = lBound + (NPoints-1)*dx; //1 point is used at each boundary
    boundaryCondition = "Periodic"; //Only boundary condition permitted at present in periodic
    gridLocations = (double*) calloc(NP, sizeof(double));
    EfieldValues = (double*) calloc(NP, sizeof(double));
    
    BfieldY = (double*) calloc(NP, sizeof(double));
    BfieldZ = (double*) calloc(NP, sizeof(double));

    chargeDensity = (double*) calloc(NP, sizeof(double));
    jValues = (double*) calloc(NP, sizeof(double));

    timestep = 0;
    phi = (double*) calloc(NP, sizeof(double));
    dt = dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
        EfieldValues[i] = 0;  
        BfieldValues[i] = 0;
        chargeDensity[i] = 0;
    }

    poissonDiagonal = (double*) calloc((NPoints-2), sizeof(double));
    poissonUpDiagonal = (double*) calloc((NPoints-3), sizeof(double));
    poissonLDiagonal = (double*) calloc(NPoints - 3, sizeof(double));
    for(int i = 0; i < NPoints-3; i++){
        poissonDiagonal[i] = -2;
        poissonUpDiagonal[i] = 1;
        poissonLDiagonal[i] = 1;
    }
    poissonDiagonal[(NPoints - 2)-1] = -2;
    
    int info = 0;
   
    MKL_INT NPl2 = NP-2;

    //dgetrf(&NPl2,&NPl2,poissonMatrix,&NPl2,ipiv,&info);
    ddttrfb(&NPl2, poissonLDiagonal,poissonDiagonal,poissonLDiagonal,&info);


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

    poissonDiagonal = (double*) calloc((NPoints-2), sizeof(double));
    poissonUpDiagonal = (double*) calloc((NPoints-3), sizeof(double));
    poissonLDiagonal = (double*) calloc(NPoints - 3, sizeof(double));
    for(int i = 0; i < NPoints-3; i++){
        poissonDiagonal[i] = -2;
        poissonUpDiagonal[i] = 1;
        poissonLDiagonal[i] = 1;
    }
    poissonDiagonal[(NPoints - 2)-1] = -2;
    
    int info = 0;
   
    MKL_INT NPl2 = NP-2;

    //dgetrf(&NPl2,&NPl2,poissonMatrix,&NPl2,ipiv,&info);
    ddttrfb(&NPl2, poissonLDiagonal,poissonDiagonal,poissonLDiagonal,&info);

}

void Grid<1>::vertexInfoTraverse(){
    std::ofstream ofile;
    std::string filename = "dumps/gridQuants" + std::to_string(dump) + ".txt";
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
    std::string filename = "dumps/particleProps" + std::to_string(dump) + ".txt";
    ofile.open(filename);
    ofile << "posits " << "vels " << std::endl;
    for(int i = 0; i < NParticles; i++){
        ofile << pList[i].position << " ";
        ofile << pList[i].velocity << std::endl;
    }
    ofile.close();
}

void Grid<1>::poissonSolver(){
    const char trans = 'N';
    double* density = (double*) malloc(sizeof(double) * (NP-2));
    for(int i = 0; i < NP-2; i++){
        density[i] = -1*chargeDensity[i+1] * dx * dx;
        //std::cout << "Density vector at " << i << " is: " << density[i] << std::endl;
    }
   
    int info = 0;
    MKL_INT nrhs = 1;
    MKL_INT NPl2 = NP-2;

    
    ddttrsb(&trans,&NPl2,&nrhs,poissonLDiagonal,poissonDiagonal,poissonUpDiagonal,density,&NPl2,&info);
    
    
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
    double rdx = 1.0/dx;
    
#pragma omp parallel for default(none) shared(rdx, pList, NParticles), schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        int lIndex = pList[i].lIndex;
        int rIndex = pList[i].rIndex;
        
        double effectiveEfield = EfieldValues[lIndex] * ((dx*rIndex)- pList[i].position) * rdx + 
                                EfieldValues[rIndex] * (pList[i].position - (dx*lIndex)) * rdx;
        pList[i].velocity -= pList[i].qmRatio * effectiveEfield * (dt/2.0);
    }
}

void Grid<1>::updateVelocities(Particle<1>* pList, int NParticles){
   
    double rdx = 1.0/dx;
    
#pragma omp parallel for default(none) shared(rdx, pList, NParticles), schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        int lIndex = pList[i].lIndex;
        int rIndex = pList[i].rIndex;
        double effectiveEfield = EfieldValues[lIndex] * ((dx*rIndex)- pList[i].position) * rdx + 
                                EfieldValues[rIndex] * (pList[i].position - (dx*lIndex)) * rdx;
        pList[i].velocity += pList[i].qmRatio * effectiveEfield * (dt);
    }
        //rotate velocities
        
}

void Grid<1>::moveParticles(Particle<1>* pList, int NParticles){
    double rdx = 1.0/dx;
#pragma omp parallel for default(none) shared(NParticles, pList, rdx) schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        pList[i].position += pList[i].velocity * dt;
        if(pList[i].position > rBound) 
            pList[i].position = lBound + (pList[i].position-rBound);
        if(pList[i].position < lBound)
            pList[i].position = rBound - (lBound - pList[i].position);
        
        pList[i].lIndex = ((double) pList[i].position-lBound) * rdx;
        pList[i].rIndex = (pList[i].lIndex + 1);

        if(pList[i].rIndex > NP-1){
            pList[i].rIndex = 0;
        }
    }
    timestep++;
    time += dt;
}

void Grid<1>::Initialize(Particle<1>* pList, const int NParticles){
    
    particleInitUniformProtonElectronPairs(pList, NParticles);
    
    particleInCell(pList, NParticles);
    
    auto tStart = std::chrono::high_resolution_clock::now();
    poissonSolver();
    auto tEnd = std::chrono::high_resolution_clock::now();
    auto tTime =  std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
    std::cout << "Poisson Solver took: " << tTime.count() << std::endl;

    vertexInfoTraverse();

    getInitialVelocities(pList, NParticles);

    particleInfoTraverse(pList, NParticles);

    dump++;

}


void Grid<1>::IntegrationLoop(Particle<1>* pList, const int NParticles){
        auto loopStart = std::chrono::high_resolution_clock::now();
        updateVelocities(pList,NParticles);
        auto velUpdate = std::chrono::high_resolution_clock::now();
        moveParticles(pList, NParticles);
        auto moveUpdate = std::chrono::high_resolution_clock::now();
        particleInCell(pList, NParticles);
        auto pInCell = std::chrono::high_resolution_clock::now();
        //auto tStart = std::chrono::high_resolution_clock::now();
        poissonSolver();
        //auto tEnd = std::chrono::high_resolution_clock::now();
        //auto tTime =  std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
        //std::cout << "Poisson Solver took: " << tTime.count() << std::endl;


        if(timestep % timeStepRate == 0){
            vertexInfoTraverse();
            particleInfoTraverse(pList, NParticles);
            dump++;
        }
        auto loopEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Int loop: " << 
                std::chrono::duration_cast<std::chrono::duration<double>>(loopEnd-loopStart).count() << "\t";
        std::cout << "velUp: " << std::chrono::duration_cast<std::chrono::duration<double>>(velUpdate-loopStart).count()<< "\t";
        std::cout << "posUp: " << std::chrono::duration_cast<std::chrono::duration<double>>(moveUpdate-velUpdate).count() << "\t";
        std::cout << "pInCell: " << std::chrono::duration_cast<std::chrono::duration<double>>(pInCell - moveUpdate).count() << "\t";
        std::cout << "poisson: " << std::chrono::duration_cast<std::chrono::duration<double>>(loopEnd - pInCell).count() << "\t";
        std::cout << "time: " << time << std::endl;
}