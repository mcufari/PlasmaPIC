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
    const int Nparticles = 2;
    const double tMax = 0.1;
    Grid<1> g(11, 0.1, 0);
    
    Particle<1>* pList = (Particle<1>*) malloc(sizeof(Particle<1>) * Nparticles);
    
    g.Initialize(pList, Nparticles);
    
    std::cout << "Time stepping at: " << g.dt << std::endl;
    double t = 0;
    while(t < tMax){
        g.IntegrationLoop(pList, Nparticles);
        t+= g.dt;
    }
    
    
}