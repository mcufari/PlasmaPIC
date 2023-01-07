#include "grid/Grid.hpp"
#include "particles/Particle.hpp"
#include "particles/particleInCell.hpp"

int main(){
    //Read in setup routine and dynamically select routine
    //Read in initial conditions
    //Read in boundary conditions
    //Initialize problem and run tests
    //Output init
    //Enter simulation loop
    //Exit simulation loop
    //Dump info
    //Dump closing remarks
    //End
    const int Nparticles = 10;
    Grid<1> g(10, 5, 0);
    g.vertexInfoTraverse();
    Particle<1>* pList = (Particle<1>*) malloc(sizeof(Particle<1>) * Nparticles);
    particleInitRandomStaticProton(pList, Nparticles, g);
    particleInCell(pList, Nparticles, g);
    g.vertexInfoTraverse();
}