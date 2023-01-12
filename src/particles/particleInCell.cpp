#include "../grid/Grid.hpp"
#include <math.h>
#include <string.h>
#include "omp.h"
void Grid<1>::particleInCell(const Particle<1>* pList, const int NParts) {
    memset(&chargeDensity[0], 0, NP*sizeof(double));
    
    double rdx = 1.0/dx;
    for(int i = 0; i < NParts; i++){
        int lIndex = pList[i].lIndex;
        int rIndex = pList[i].rIndex;
        pList[i].cloudInCell(&chargeDensity[lIndex],&chargeDensity[rIndex],lIndex*dx,rIndex*dx, rdx);
    }
}

void Grid<1>::currentInCell(const Particle<1>* pList, const int NParts, double* currentList){
    double rdx = 1.0/dx;
    for(int i = 0; i < NParts; i++){
        int lhs = pList[i].lIndex ;
        int rhs = pList[i].rIndex;
        std::cout << "rhs: " << rhs << std::endl;
        std::cout << "lhs: " << lhs << std::endl;
        pList[i].currentCloudInCell(&currentList[lhs], &currentList[rhs], lhs*dx, rhs*dx, rdx);
    }
}

void Grid<1>::particleInitRandomStaticProton(Particle<1>* pList, const int NParts) const {
    for(int i = 0; i < NParts; i++){
        pList[i] = Particle<1>(1,0.15,0,0.01);
   
        pList[i].velocity = 0;
      
        //pList[i].position = ((double) rand())/(double) RAND_MAX * (rBound - lBound) + lBound;
        //std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
}

void Grid<1>::particleInitStaticProtonPair(Particle<1>* pList, const int NParts){
    if(NParts % 2 == 1){
        std::cout << "Error, nparts must be an even number" << std::endl;
        exit(1);
    }
    for(int i = 0; i < NParts/2; i++){
        pList[i] = Particle<1>(-1,(double) rand() / RAND_MAX, 0, 0.01);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
        //std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
    for(int i = NParts/2; i < NParts; i++){
        pList[i] = Particle<1>(1, (double) rand() / RAND_MAX, 0 ,0.01);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
        //std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
    

}

void Grid<1>::particleInitUniformProtonElectronPairs(Particle<1>* pList, const int NParts){
    for(int i = 0; i < NParts; i+=2){
        pList[i] = Particle<1>(-1, (double) (i+1) /(NParts+1), 0, 0.01);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
    }
    for(int i = 1; i < NParts; i += 2){
        pList[i] = Particle<1>(1, (double)(i+1) /(NParts+1), 0, 1.0);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
    }
}

void Grid<1>::particleInitSinusoidElectrons(Particle<1>* pList, const int NParts){
    for(int i = 0; i < NParts/2; i++){
        pList[i] = Particle<1>(-1, 0, 0, 0.01);
        double pos = 1.0/(2.0*M_PI) * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX );
        //std::cout << "placing particle at: " << pos << std::endl;
        pList[i].position= pos;
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;

    }
    for(int i = NParts/2; i < NParts; i++){
        pList[i] = Particle<1>(1, 0, 0, 1.0);
        double pos = 1.0/2.0 + 1.0/(2.0*M_PI) * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX );
        //std::cout << "placing particle at: " << pos << std::endl;
        pList[i].position= pos;
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
    }
    
}