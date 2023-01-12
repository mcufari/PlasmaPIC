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
    
    EFieldY = (double*) calloc(NP, sizeof(double));
    EFieldZ = (double*) calloc(NP, sizeof(double));
    fzDownValues = (double*) calloc(NP, sizeof(double));
    fzUpValues = (double*) calloc(NP, sizeof(double));

   
    BFieldY = (double*) calloc(NP, sizeof(double));
    BFieldZ = (double*) calloc(NP, sizeof(double));

    chargeDensity = (double*) calloc(NP, sizeof(double));
    jyUpValues = (double*) calloc(NP, sizeof(double));
    jyDownValues = (double*) calloc(NP, sizeof(double));
    fzUpValues = (double*) calloc(NP, sizeof(double));
    fzDownValues = (double*) calloc(NP, sizeof(double));

    timestep = 0;
    phi = (double*) calloc(NP, sizeof(double));
    dt = dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
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



/*Grid<1>::Grid(int NPoints, double dx, double lBoundary, std::string boundCondition) : 
NP(NPoints), 
dx(dx), 
lBound(lBoundary), 
boundaryCondition(boundCondition) 
{
    rBound = lBound + (NPoints - 1) * dx;
    gridLocations = (double*) calloc(NP, sizeof(double));
    EfieldValues = (double*) calloc(NP, sizeof(double));
    BfieldValues = (double*) calloc(NP, sizeof(double));
    
    EFieldY = (double*) calloc(NP, sizeof(double));
    EFieldZ = (double*) calloc(NP, sizeof(double));
    fzDownValues = (double*) calloc(NP, sizeof(double));
    fzUpValues = (double*) calloc(NP, sizeof(double));

    BFieldY = (double*) calloc(NP, sizeof(double));
    BFieldZ = (double*) calloc(NP, sizeof(double));

    chargeDensity = (double*) calloc(NP, sizeof(double));
    jyUpValues = (double*) calloc(NP, sizeof(double));
    jyDownValues = (double*) calloc(NP, sizeof(double));
    
    timestep = 0;
    phi = (double*) calloc(NP, sizeof(double));
    dt = dx;
    for(int i = 0; i < NPoints; i++){
        gridLocations[i] = (lBound + i*dx);
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

}*/

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
    ofile << "posits " << "vels " <<  "lIndex " << "rIndex " << std::endl;
    for(int i = 0; i < NParticles; i++){
        ofile << pList[i].position << " ";
        ofile << pList[i].velocity << " ";
        ofile << pList[i].lIndex << " ";
        ofile << pList[i].rIndex << std::endl;
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

void Grid<1>::EYBZSolver(){
    double* fzUpValuesNew = (double*) calloc(NP,sizeof(double));
    double* fzDownValuesNew = (double*) calloc(NP, sizeof(double));
    std::cout << "allocated fzDown new and fzupnew" << std::endl;

    for(int i = 1; i < NP-1; i++){
        fzUpValuesNew[i] = fzUpValues[i-1] - dt/4.0 * (jyDownValues[i-1] + jyUpValues[i]);
        fzDownValuesNew[i] = fzDownValues[i+1] - dt/4.0 * (jyDownValues[i+1] + jyUpValues[i]);
    }
    std::cout << "finished first iter" << std::endl;

    fzUpValuesNew[0] = fzUpValues[NP-1] - dt/4.0 * (jyDownValues[NP-1] + jyUpValues[0]);
    fzDownValuesNew[0] = fzDownValues[1] - dt/4.0 * (jyDownValues[1] + jyDownValues[0]);

    fzUpValuesNew[NP-1] = fzUpValues[NP-2] - dt/4.0 * (jyDownValues[NP-2] + jyUpValues[NP-1]);
    fzDownValuesNew[NP-1] = fzDownValues[0] - dt/4.0 * (jyDownValues[0] + jyUpValues[NP-1]);
    std::cout << "finished boundary" << std::endl;

    for(int i = 0; i < NP; i++){
        std::cout << i << std::endl;
        std::cout << "The value of EFieldY[i]: " << EFieldY[i] << std::endl;
        std::cout << "The value of BFieldZ[i]: " << BFieldZ[i] << std::endl;
        std::cout << "the value of fzUpValuesNew[i]: " << fzUpValuesNew[i] << std::endl;
        std::cout << "the value of fzDownValuesNew[i]: " << fzDownValuesNew[i] << std::endl;
        EFieldY[i] = fzUpValuesNew[i] + fzDownValuesNew[i];
        std::cout << "Efield update complete" << std::endl;
        BFieldZ[i] = fzUpValuesNew[i] - fzDownValuesNew[i];
        std::cout << "BField update complete" << std::endl;
    }
    std::cout << "Finished EField BField update" << std::endl;
    memcpy(fzUpValues, fzUpValuesNew, sizeof(double)*NP);
    memcpy(fzDownValues, fzDownValuesNew, sizeof(double)*NP);
    std::cout << "finished memcopy" << std::endl;
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
    
//#pragma omp parallel for default(none) shared(rdx, pList, NParticles), schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        int lIndex = pList[i].lIndex;
        int rIndex = pList[i].rIndex;
        double pos = pList[i].position;
        double vx = pList[i].velocity;
        double vy = pList[i].vely;
        double vz = pList[i].velz;

        double qPrime = pList[i].qmRatio * dt/2.0;

        double lweight = (dx*rIndex - pos)*rdx;
        double rweight = (pos - dx*lIndex)*rdx;
        double Ex = EfieldValues[lIndex] * lweight  + 
                                EfieldValues[rIndex] * rweight ;
        double Ey = EFieldY[lIndex] * lweight  + 
                                EFieldY[rIndex] * rweight;
        double Ez = EFieldZ[lIndex] * lweight  + 
                                EFieldZ[rIndex] * rweight;
        
        double Bx = 0;
        double By = BFieldY[lIndex] * lweight  + BFieldY[rIndex]*rweight;
        double Bz = BFieldZ[lIndex] * lweight + BFieldZ[rIndex] * rweight;

        double ux = vx + qPrime * Ex;
        double uy = vy + qPrime * Ey;
        double uz = vz + qPrime * Ez;

        double hx = qPrime * Bx;
        double hy = qPrime * By;
        double hz = qPrime * Bz;

        double h2 = hx*hx + hy*hy + hz*hz;
        double sx = 2.0 * hx / (1.0 + h2);
        double sy = 2.0 * hy / (1.0 + h2);
        double sz = 2.0 * sz / (1.0 + h2);

        double uPrimeX = ux - hy*sy*ux - hz*sz*ux + hx*sy*uy + sz*uy - sy*uz + hx*sz*uz;
        double uPrimeY = hy*sx*ux - sz*ux + uy - hx*sx*uy - hz*sz*uy + sx*uz + hy*sz*uz;
        double uPrimeZ = hz*sx*ux + sy*ux - sx*uy + hz*sy*uy + uz - hx*sx*uz - hy*sy*uz;
        pList[i].velocity = uPrimeX + qPrime * Ex;
        pList[i].vely = uPrimeY + qPrime*Ey;
        pList[i].velz = uPrimeZ + qPrime*Ez;

    }
        //rotate velocities
    
}

void Grid<1>::moveParticles(Particle<1>* pList, int NParticles){
    double rdx = 1.0/dx;
//#pragma omp parallel for default(none) shared(NParticles, pList, rdx) schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        pList[i].oldPosition = pList[i].position;
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
    std::cout << "initalizing" << std::endl;
    particleInitUniformProtonElectronPairs(pList, NParticles);
    std::cout << "initalized particle locations" << std::endl;
    particleInCell(pList, NParticles);
    std::cout << "performed particle in cell";
    auto tStart = std::chrono::high_resolution_clock::now();
    poissonSolver();
    std::cout << "performed poisson solver" << std::endl;
    auto tEnd = std::chrono::high_resolution_clock::now();
    auto tTime =  std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
    std::cout << "Poisson Solver took: " << tTime.count() << std::endl;
    
    vertexInfoTraverse();
    std::cout << "performed vertex traverse" << std::endl;
    getInitialVelocities(pList, NParticles);
    std::cout << "got initial vs" << std::endl;
    particleInfoTraverse(pList, NParticles);
    std::cout << "got particle info" << std::endl;
    dump++;

}


void Grid<1>::IntegrationLoop(Particle<1>* pList, const int NParticles){
        auto loopStart = std::chrono::high_resolution_clock::now();
        //Get v(t+1/2*dt);
        updateVelocities(pList,NParticles);
        std::cout << "updated vels" << std::endl;
        currentInCell(pList, NParticles, jyDownValues);
        std::cout << "updated current in cells" << std::endl;
        auto velUpdate = std::chrono::high_resolution_clock::now();
        moveParticles(pList, NParticles);
        std::cout << "particles moved" << std::endl;
        //get x(t+dt)
        auto moveUpdate = std::chrono::high_resolution_clock::now();
        currentInCell(pList, NParticles, jyUpValues);
        std::cout << "current in cell calculated" << std::endl;
        //update currents at J(t+dt)
        particleInCell(pList, NParticles);
        std::cout << "particle in cell calculated" << std::endl;
     

        auto pInCell = std::chrono::high_resolution_clock::now();
        //auto tStart = std::chrono::high_resolution_clock::now();
        poissonSolver();

        std::cout << "particle in cell performed" << std::endl;
        EYBZSolver();
        std::cout << "Efield solver performed" << std::endl;
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