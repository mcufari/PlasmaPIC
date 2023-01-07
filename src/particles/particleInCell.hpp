#ifndef PARTICLE_IN_CELL_SUBROUTINES
#define PARTICLE_IN_CELL_SUBROUTINES

#include "Particle.hpp"
#include "../grid/Grid.hpp"
void particleInCell(const Particle<1>* pList, const int NParts, Grid<1> grid); //Traverse the pList and update the corresponding gridVertices 
void particleInitRandomStaticProton(Particle<1>* pList, const int NParts, const Grid<1> grid); //Place particles at random locations on the grid
#endif