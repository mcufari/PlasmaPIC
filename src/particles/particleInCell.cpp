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
        pList[i] = Particle<1>(1,.15,0,1);
   
        pList[i].velocity = 0;
      
        //pList[i].position = ((double) rand())/(double) RAND_MAX * (grid.rBound - grid.lBound) + grid.lBound;
        std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
}

void Grid<1>::particleInitStaticProtonPair(Particle<1>* pList, const int NParts){
    pList[0] = Particle<1>(1, 0.125, 0, 1);
    pList[1] = Particle<1>(1, 0.175, 0, 1);
}
