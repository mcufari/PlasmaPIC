#include "../grid/Grid.hpp"
#include <math.h>
#include <string.h>
#include "omp.h"
void Grid<1>::particleInCell(double * pList, int* iList, const int NParts) {
    memset(&chargeDensity[0], 0, NP*sizeof(double));
    
    double rdx = 1.0/dx;
    for(int i = 0; i < NParts; i++){
        int lIndex = iList[i];
        int rIndex = iList[NParts + i];
        chargeDensity[lIndex] = pList[i] * ((rIndex*dx) - pList[3*NParts + i]) * rdx;
    
        chargeDensity[rIndex] = pList[i] * (pList[3*NParts + i] - (lIndex*dx)) * rdx;

    }
}

void Grid<1>::currentInCell(double* pList, int* iList, const int NParts, double* j){
    double rdx = 1.0/dx;
    for(int i = 0; i < NParts; i++){
        int lhs = iList[i] ;
        int rhs = iList[NParts + i];
        
        j[lhs] = pList[i] * pList[5*NParts + i] * ((rhs*dx) - pList[3*NParts + i]) * rdx;
    
        j[rhs] = pList[i] * pList[5*NParts + i] * (pList[3*NParts + i] - (lhs*dx)) * rdx;

    }
}

void Grid<1>::particleInitRandomStaticProton(Particle<1>* pList, const int NParts) const {
    for(int i = 0; i < NParts; i++){
        pList[i] = Particle<1>(1,0.15,0,0.01);
   
        pList[i].velocity = 0;
      
       
    }
}

void Grid<1>::particleInitStaticProtonPair(Particle<1>* pList, const int NParts){
    if(NParts % 2 == 1){
        std::cout << "Error, nparts must be an even number" << std::endl;
        exit(1);
    }
    for(int i = 0; i < NParts/2; i++){
        pList[i] = Particle<1>(-1, 0.1756 , 0, 0.01);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
        
    }
    for(int i = NParts/2; i < NParts; i++){
        pList[i] = Particle<1>(1, 0.154, 0 ,0.01);
        pList[i].lIndex = pList[i].position / dx;
        pList[i].rIndex = (pList[i].lIndex + 1) % NP;
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

void Grid<1>::particleInitSinusoidElectrons(double* pList, int* indexAr, const int NParts){
    for(int i = 0; i < NParts/2; i++){
        pList[i] = -1.0; // charge
        pList[NParts + i] = 1.0; // mass 
        pList[2*NParts + i] = -1.0; // q/M
        double pos = 1.0/M_PI * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX ); 
        pList[3*NParts + i] = 1.0/M_PI * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX ); // x
        pList[4*NParts + i] = 0.0; // vx
        pList[5*NParts + i] = 0.0; //vy
        pList[6*NParts + i] = 0.0; // vz
        indexAr[i] =(int) pos / dx;
        indexAr[NParts + i] =(int) (pos/dx + 1) % NP;

    }
    for(int i = NParts/2; i < NParts; i++){
        pList[i] = 1.0; // charge
        pList[NParts + i] = 1.0; // mass 
        pList[2*NParts + i] = 1.0; // q/M
        double pos = 1.0/M_PI * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX ); 
        pList[3*NParts + i] = 1.0/M_PI * acos(1.0 - 2.0 * (double) rand() / (double) RAND_MAX ); // x
        pList[4*NParts + i] = 0.0; // vx
        pList[5*NParts + i] = 0.0; //vy
        pList[6*NParts + i] = 0.0; // vz
        indexAr[i] = (int) pos / dx;
        indexAr[NParts + i] = (int) (pos/dx + 1) % NP;

    }
    
}