#include "Particle.hpp"

Particle<1>::Particle(double ch, double pos, double vel): 
    charge(ch), 
    position(pos),
    velocity(vel)
{}

void Particle<1>::cloudInCell(double* lhCharge, double* rhCharge, double* lhs, double* rhs, const double dx) const {
   *lhCharge += this->charge * ((*rhs) - this->position)/dx;
    *rhCharge += this->charge * (this->position - (*lhs))/dx;
}