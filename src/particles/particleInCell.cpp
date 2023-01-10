#include "../grid/Grid.hpp"
#include <math.h>

void Grid<1>::particleInCell(const Particle<1>* pList, const int NParts) {
    for(int i = 0; i < NP; i++){
        chargeDensity[i] = 0;
    }

    for(int i = 0; i < NParts; i++){
        int lIndex = floor((pList[i].position - lBound)/dx);
        int rIndex = (lIndex + 1) % NP;
        pList[i].cloudInCell(&chargeDensity[lIndex],&chargeDensity[rIndex],&gridLocations[lIndex],&gridLocations[rIndex],dx);
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
        pList[i] = Particle<1>(-1,lBound + dx + ((rBound-dx) - (lBound + dx))*((double)rand())/RAND_MAX, 0, 0.01);
        //std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
    for(int i = NParts/2; i < NParts; i++){
        pList[i] = Particle<1>(1, lBound + dx + ((rBound-dx) - (lBound + dx))*((double)rand())/RAND_MAX, 0 ,0.01);
        //std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
    

}

void Grid<1>::particleInitUniformProtonElectronPairs(Particle<1>* pList, const int NParts){
    for(int i = 0; i < NParts; i+=2){
        pList[i] = Particle<1>(-1, lBound + dx + ((rBound-dx) - (lBound + dx)) * i / NParts, 0, 0.01);
    }
    for(int i = 1; i < NParts; i += 2){
        pList[i] = Particle<1>(1, lBound + dx + ((rBound-dx) - (lBound + dx)) * i / NParts, 0, 0.01);
    }
}
