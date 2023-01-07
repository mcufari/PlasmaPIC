
#ifndef PARTICLE_HEADER
#define PARTICLE_HEADER
#include "../grid/GridVertex.hpp"

template<int dim = 0>
class Particle{
public:
    Particle();
};

template<> class Particle<1>{
public:
    Particle(double ch, double pos, double vel);
    void cloudInCell(GridVertex<1>* gvl, GridVertex<1>* gvr, const double dx) const; 
    double charge;
    double position;
    double velocity;
};



#endif 