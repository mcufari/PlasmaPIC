
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
    Particle(const Particle<1>& prhs);
    void cloudInCell(double* lhCharge, double* rhCharge, double lhs, double rhs, const double rdx) const;
    void currentCloudInCell(double* lhCurrent, double* rhCurrent, double lhs, double rhs, const double rdx) const;
    double charge;
    double mass;
    double qmRatio;
    
    double position;
    
    double oldPosition;
    double velocity;
    

    int lIndex;
    int rIndex;

    double vely;
    double velz;

};



#endif 