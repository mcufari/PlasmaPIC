#include "Particle.hpp"

Particle<1>::Particle(double ch, double pos, double vel): 
    charge(ch), 
    position(pos),
    velocity(vel)
{}

void Particle<1>::cloudInCell(GridVertex<1>* gvl, GridVertex<1>* gvr, const double dx) const {
    gvl->effectiveCharge += this->charge * (gvr->position - this->position)/dx;
    gvr->effectiveCharge += this->charge * (this->position - gvl->position)/dx;
}