
#ifndef PARTICLE_HEADER
#define PARTICLE_HEADER


template<int dim = 0>
class Particle{
public:
    Particle();
};

template<> class Particle<1>{
public:
    Particle(double ch, double pos, double vel, double mass);
    void cloudInCell(double* lhCharge, double* rhCharge, double* lhs, double* rhs, const double dx) const;
    double charge;
    double position;
    double velocity;
    double mass;

};



#endif 