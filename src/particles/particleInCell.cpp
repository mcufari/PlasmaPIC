#include "particleInCell.hpp"
#include <math.h>

void particleInCell(const Particle<1>* pList, const int NParts, Grid<1> grid){
    
    for(int i = 0; i < NParts; i++){
        int lIndex = floor((pList[i].position - grid.lBound)/grid.dx);
        int rIndex = (lIndex + 1) % grid.NP;
        pList[i].cloudInCell(&grid.chargeDensity[lIndex],&grid.chargeDensity[rIndex],&grid.gridLocations[lIndex],&grid.gridLocations[rIndex],grid.dx);
    }

}

void particleInitRandomStaticProton(Particle<1>* pList, const int NParts, const Grid<1> grid){
    for(int i = 0; i < NParts; i++){
        pList[i].charge = 1;
        pList[i].velocity = 0;
        pList[i].position = ((double) rand())/(double) RAND_MAX * (grid.rBound - grid.lBound) + grid.lBound;
        std::cout << "Placing particle at: " << pList[i].position << std::endl;
    }
}
