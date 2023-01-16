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
    ofile << "gridPosits " << "chargeDensity " << "Ex Ey Ez"  << " phi " << "By " << "Bz" << std::endl;
    for(int i = 0; i < NP; i++){
        ofile << gridLocations[i] << " ";
        ofile << chargeDensity[i] << " ";
        ofile << EfieldValues[i] << " " << EFieldY[i] << " " << EFieldZ[i] << " ";
        ofile << phi[i] << " ";
        ofile << BFieldY[i] << " " << BFieldZ[i] << std::endl;

    }
    ofile.close();
}

void Grid<1>::particleInfoTraverse(double* pList, int* iList, int NParticles){
    std::ofstream ofile;
    std::string filename = "dumps/particleProps" + std::to_string(dump) + ".txt";
    ofile.open(filename);
    ofile << "x " << "vx vy vz" <<  "lIndex " << "rIndex " << std::endl;
    for(int i = 0; i < NParticles; i++){
        ofile << pList[3*NParticles + i] << " ";
        ofile << pList[4*NParticles+i] << " " << pList[5*NParticles + i] << " " << pList[NParticles + i]<< " ";
        ofile << iList[i] << " " << iList[NParticles + i] << std::endl;
    }
    ofile.close();
}

void Grid<1>::poissonSolver(){
    const char trans = 'N';
    double* density = (double*) malloc(sizeof(double) * (NP-2));
    for(int i = 0; i < NP-2; i++){
        density[i] = -1*chargeDensity[i+1] * dx * dx;
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

    EfieldValues[0] = (-1.0*phi[1])/(dx);
    EfieldValues[NP-1] = (phi[NP-2])/(dx);
    for(int i = 1; i < NP-1; i++){
        EfieldValues[i] = (phi[i-1] - phi[i+1])/(2.0*dx);
    }
    
   
    free(density);
}

void Grid<1>::EYBZSolver(){
    double* fzUpValuesNew = (double*) calloc(NP,sizeof(double));
    double* fzDownValuesNew = (double*) calloc(NP, sizeof(double));

    for(int i = 1; i < NP-1; i++){
        fzUpValuesNew[i] = fzUpValues[i-1] - dt/4.0 * (jyDownValues[i-1] + jyUpValues[i]);
        fzDownValuesNew[i] = fzDownValues[i+1] - dt/4.0 * (jyDownValues[i+1] + jyUpValues[i]);
    }

    fzUpValuesNew[0] = fzUpValues[NP-1] - dt/4.0 * (jyDownValues[NP-1] + jyUpValues[0]);
    fzDownValuesNew[0] = fzDownValues[1] - dt/4.0 * (jyDownValues[1] + jyDownValues[0]);

    fzUpValuesNew[NP-1] = fzUpValues[NP-2] - dt/4.0 * (jyDownValues[NP-2] + jyUpValues[NP-1]);
    fzDownValuesNew[NP-1] = fzDownValues[0] - dt/4.0 * (jyDownValues[0] + jyUpValues[NP-1]);

    for(int i = 0; i < NP; i++){
       
        EFieldY[i] = fzUpValuesNew[i] + fzDownValuesNew[i];
        BFieldZ[i] = fzUpValuesNew[i] - fzDownValuesNew[i];
    
    }
    
    memcpy(fzUpValues, fzUpValuesNew, sizeof(double)*NP);
    memcpy(fzDownValues, fzDownValuesNew, sizeof(double)*NP);
   
}

void Grid<1>::getInitialVelocities(double* pList, int* iList, int NParticles){
    double rdx = 1.0/dx;
    
#pragma omp parallel for default(none) shared(rdx, pList, iList, NParticles), schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        int lIndex = iList[i];
        int rIndex = iList[NParticles + i];
        
        double effectiveEfield = EfieldValues[lIndex] * ((dx*rIndex)- pList[3*NParticles + i]) * rdx + 
                                EfieldValues[rIndex] * (pList[3*NParticles + i] - (dx*lIndex)) * rdx;
        pList[4*NParticles + i] -= pList[2 * NParticles + i] * effectiveEfield * (dt/2.0);
    }
}

void Grid<1>::updateVelocities(double* pList, int* iList, int NParticles){
    double rdx = 1.0/dx;
    size_t i = 0;
    
    double dtAr[4] = {dt/2.0, dt/2.0, dt/2.0, dt/2.0};
    __m256d dt2 = _mm256_load_pd(&dtAr[0]);
    double lw;
    double rw;
    double efx[4];
    double efy[4];
    double efz[4];
    double hfx[4];
    double hfy[4];
    double hfz[4];
    double hfsq[4];

    for(; i < (NParticles & (~0x3)); i+=4){

        __m256d pxd = _mm256_load_pd(&pList[3*NParticles + i]);
        __m256d vxa = _mm256_load_pd(&pList[4*NParticles + i]);
        __m256d vya = _mm256_load_pd(&pList[5*NParticles + i]);
        __m256d vza = _mm256_load_pd(&pList[6*NParticles + i]);
        __m256d qPa = _mm256_load_pd(&pList[7*NParticles + i]);
        for(int j = 0; j < 4; j++){
            lw = iList[i+j]*dx - pList[3*NParticles + i + j];
            rw = pList[3*NParticles +i + j] - iList[NParticles + i+j]*dx;
            efx[j] = EfieldValues[iList[i+j]] * lw + EfieldValues[iList[NParticles+i+j]]*rw;
            efy[j] = EFieldY[iList[i+j]] * lw + EFieldY[iList[NParticles+i+j]]*rw;
            efz[j] = EFieldZ[iList[i+j]] * lw + EFieldZ[iList[NParticles+i+j]]*rw;
            
            hfy[j] = BFieldY[iList[i+j]] * lw + BFieldY[iList[NParticles+i+j]]*rw;
            hfz[j] = BFieldZ[iList[i+j]] * lw + BFieldZ[iList[NParticles+i+j]]*rw;
            hfsq[j] = 2.0 / (1.0 + hfy[j]*hfy[j] + hfz[j]*hfz[j]);
        }
       
        qPa = _mm256_mul_pd(qPa, dt2);
        __m256d edx = _mm256_load_pd(&efx[0]);
        __m256d edy = _mm256_load_pd(&efy[0]);
        __m256d edz = _mm256_load_pd(&efz[0]);
        
        __m256d hdy = _mm256_load_pd(&hfy[0]);
        __m256d hdz = _mm256_load_pd(&hdz[0]);

        edx = _mm256_mul_pd(edx, qPa);
        edy = _mm256_mul_pd(edy, qPa);
        edz = _mm256_mul_pd(edz, qPa);
        hdy = _mm256_mul_pd(hdy, qPa);
        hdz = _mm256_mul_pd(hdz, qPa);

        __m256d udx = _mm256_add_pd(vxa, edx);
        __m256d udy = _mm256_add_pd(vya, edy);
        __m256d udz = _mm256_add_pd(vza, edz);
        __m256d sdy = _mm256_load_pd(&hfsq[0]);
        sdy = _mm256_mul_pd(sdy, hdy);
        __m256d sdz = _mm256_load_pd(&hfsq[0]);
        sdz = _mm256_mul_pd(sdz, hdz);

        __m256d temp1 = _mm256_sub_pd(udx, edx);
        __m256d temp2 = _mm256_mul_pd(hdy, sdy);
        temp2 = _mm256_mul_pd(temp2, udx);
        __m256d temp3 = _mm256_mul_pd(hdz, sdz);
        temp3 = _mm256_mul_pd(temp3, udx);
        __m256d temp4 = _mm256_mul_pd(sdz, udy);
        __m256d temp5 = _mm256_mul_pd(sdy, udz);
        temp1 = _mm256_sub_pd(temp1, temp2);
        temp1 = _mm256_sub_pd(temp1, temp3);
        temp1 = _mm256_add_pd(temp1, temp4);
        temp1 = _mm256_sub_pd(temp1, temp5);

        _mm256_stream_pd(&pList[4*NParticles+i],temp1);
        temp1 = _mm256_add_pd(udy, edy);
        temp2 = _mm256_mul_pd(sdz, udx);
        temp3 = _mm256_mul_pd(hdz, sdz);
        temp3 = _mm256_mul_pd(temp3, udy);

        temp4 = _mm256_mul_pd(hdy, sdz);
        temp4 = _mm256_mul_pd(temp4, udz);

        temp1 = _mm256_sub_pd(temp1, temp2);
        temp1 = _mm256_sub_pd(temp1, temp3);
        temp1 = _mm256_add_pd(temp1, temp4);

        _mm256_stream_pd(&pList[5*NParticles + i], temp1);

        temp1 = _mm256_add_pd(udz, edz);
        temp2 = _mm256_mul_pd(sdy, udx);
        temp3 = _mm256_mul_pd(hdz, sdy);
        temp3 = _mm256_mul_pd(temp3, udy);
        temp4 = _mm256_mul_pd(hdy, sdy);
        temp4 = _mm256_mul_pd(temp4, udz);
        temp1 = _mm256_add_pd(temp1, temp2);
        temp1 = _mm256_add_pd(temp1, temp3);
        temp1 = _mm256_sub_pd(temp1, temp4);

        _mm256_stream_pd(&pList[6*NParticles + i], temp1);
    }
    std::cout << "The value of i: " << i << std::endl;
//#pragma omp parallel for default(none) shared(rdx, pList, iList, NParticles), schedule(static) num_threads(1)
    for(; i < NParticles; i++){
        int lIndex = iList[i];
        int rIndex = iList[NParticles + i];
        double pos = pList[3*NParticles + i];
        double vx = pList[4 * NParticles + i];
        double vy = pList[5 * NParticles + i];
        double vz = pList[6 * NParticles + i];

        double qPrime = pList[2*NParticles + i] * dt/2.0;

        double lweight = (dx*rIndex - pos)*rdx;
        double rweight = (pos - dx*lIndex)*rdx;
        double Ex = (EfieldValues[lIndex] * lweight  + 
                                EfieldValues[rIndex] * rweight)*qPrime; ;
        double Ey = (EFieldY[lIndex] * lweight  + 
                                EFieldY[rIndex] * rweight)*qPrime;
        double Ez = (EFieldZ[lIndex] * lweight  + 
                                EFieldZ[rIndex] * rweight)*qPrime;
        
        double hx = 0;
        double hy = (BFieldY[lIndex] * lweight  + BFieldY[rIndex]*rweight)*qPrime;
        double hz = (BFieldZ[lIndex] * lweight + BFieldZ[rIndex] * rweight)*qPrime;

        double ux = vx +  Ex;
        double uy = vy +  Ey;
        double uz = vz +  Ez;

        double hSq = hy*hy + hz*hz;
        
        hSq = 2.0 / (1.0 + hSq);
        double sx = 0;
        double sy = hSq * hy;
        double sz = hSq * hz;

        pList[4*NParticles + i] = ux - hy*sy*ux - hz*sz*ux + hx*sy*uy + sz*uy - sy*uz + hx*sz*uz +  Ex;
        pList[5*NParticles + i] = hy*sx*ux - sz*ux + uy - hx*sx*uy - hz*sz*uy + sx*uz + hy*sz*uz + Ey;
        pList[6*NParticles + i] = hz*sx*ux + sy*ux - sx*uy + hz*sy*uy + uz - hx*sx*uz - hy*sy*uz + Ez;
    }
        //rotate velocities
    
}

void Grid<1>::moveParticles(double* pList, int* iList, int NParticles){
    double rdx = 1.0/dx;
//#pragma omp parallel for default(none) shared(NParticles, pList, rdx) schedule(static) num_threads(4)
    for(int i = 0; i < NParticles; i++){
        //pList[i].oldPosition = pList[i].position;
        pList[3*NParticles + i] += pList[4*NParticles + i] * dt;
        if(pList[3*NParticles + i] > rBound) 
            pList[3*NParticles + i] = lBound + (pList[3*NParticles + i]-rBound);
        while(pList[3*NParticles + i] < lBound){
            pList[3*NParticles + i] = rBound - (lBound - pList[3*NParticles + i]);
            }
        iList[i] =(int) ( (pList[3*NParticles+i]-lBound) * rdx );
        int rIndex = (iList[i] + 1);

        if(rIndex > NP-1){
            iList[NParticles + i] = 0;
        }
    }
    timestep++;
    time += dt;
}

void Grid<1>::Initialize(double* pList, const int NParticles, int* iList){
   
    particleInitSinusoidElectrons(pList, iList, NParticles);
    
    particleInCell(pList, iList, NParticles);
    
    auto tStart = std::chrono::high_resolution_clock::now();
    poissonSolver();
    
    auto tEnd = std::chrono::high_resolution_clock::now();
    auto tTime =  std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
    
    
    //vertexInfoTraverse();
    
    getInitialVelocities(pList, iList, NParticles);
    
    particleInfoTraverse(pList, iList, NParticles);
    
    dump++;

}


void Grid<1>::IntegrationLoop(double* pList, int* iList, const int NParticles){
        auto loopStart = std::chrono::high_resolution_clock::now();
        //Get v(t+1/2*dt);
        updateVelocities(pList, iList, NParticles);
        auto velUpdate = std::chrono::high_resolution_clock::now();
        currentInCell(pList, iList, NParticles, jyDownValues);
        auto curOne = std::chrono::high_resolution_clock::now();
        moveParticles(pList, iList, NParticles);
        //get x(t+dt)
        auto moveUpdate = std::chrono::high_resolution_clock::now();
        currentInCell(pList, iList, NParticles, jyUpValues);
        auto curTwo = std::chrono::high_resolution_clock::now();
        //update currents at J(t+dt)
        particleInCell(pList, iList, NParticles);

        auto pInCell = std::chrono::high_resolution_clock::now();
        //auto tStart = std::chrono::high_resolution_clock::now();
        poissonSolver();

        auto poissonSolved = std::chrono::high_resolution_clock::now();
       
        EYBZSolver();

        auto EYBZ = std::chrono::high_resolution_clock::now();

        //auto tEnd = std::chrono::high_resolution_clock::now();
        //auto tTime =  std::chrono::duration_cast<std::chrono::duration<double>>(tEnd - tStart);
        //std::cout << "Poisson Solver took: " << tTime.count() << std::endl;

        if(timestep % timeStepRate == 0){
            //vertexInfoTraverse();
            particleInfoTraverse(pList, iList, NParticles);
            dump++;
        }

        auto loopEnd = std::chrono::high_resolution_clock::now();
        std::cout << "Int loop: " << 
                std::chrono::duration_cast<std::chrono::duration<double>>(loopEnd-loopStart).count() << "\t";
        std::cout << "velUp: " << std::chrono::duration_cast<std::chrono::duration<double>>(velUpdate-loopStart).count()<< "\t";
        std::cout << "posUp: " << std::chrono::duration_cast<std::chrono::duration<double>>(moveUpdate-curOne).count() << "\t";
        std::cout << "pInCell: " << std::chrono::duration_cast<std::chrono::duration<double>>(pInCell - curTwo).count() << "\t";
        std::cout << "poisson: " << std::chrono::duration_cast<std::chrono::duration<double>>(poissonSolved - pInCell).count() << "\t";
        std::cout << "time: " << time << std::endl;

        std::cout << "CurOne: " << std::chrono::duration_cast<std::chrono::duration<double>>(curOne - velUpdate).count() << "\t";
        std::cout << "CurTwo: " << std::chrono::duration_cast<std::chrono::duration<double>>(curTwo - moveUpdate).count() << "\t";
        std::cout << "EYBZ: " <<   std::chrono::duration_cast<std::chrono::duration<double>>(EYBZ - poissonSolved).count() << "\t";
        std::cout << "Dump: " <<   std::chrono::duration_cast<std::chrono::duration<double>>(loopEnd - EYBZ).count() << std::endl;

}