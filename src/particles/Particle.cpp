#include "Particle.hpp"

Particle<1>::Particle(double ch, double pos, double vel, double mass): 
    charge(ch), 
    position(pos),
    velocity(vel),
    mass(mass),
    qmRatio(ch/mass)
{
    lIndex  = 0;
    rIndex =  0;
    vely = 0.0;
    velz = 0.0;
    oldPosition = 0.0;
}

void Particle<1>::cloudInCell(double* lhCharge, double* rhCharge, double lhs, double rhs, const double rdx) const {
   
   *lhCharge += charge * ((rhs) - position) * rdx;
    *rhCharge += charge * (position - (lhs)) * rdx;
}

void Particle<1>::currentCloudInCell(double* lhCurrent, double* rhCurrent, double lhs, double rhs, const double rdx) const {
   
   *lhCurrent += charge * vely * ((rhs) - position) * rdx;
    *rhCurrent += charge * vely * (position - (lhs)) * rdx;
}

Particle<1>::Particle(const Particle<1>& prhs){
    mass = prhs.mass;
    charge = prhs.charge;
    velocity = prhs.velocity;
    position = prhs.position;
    qmRatio = charge/mass;
    vely = prhs.vely;
    velz = prhs.velz;
    oldPosition = prhs.oldPosition;
}