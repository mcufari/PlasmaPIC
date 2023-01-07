#ifndef GRIDVERTEX_HEADER
#define GRIDVERTEX_HEADER
//Declaration of GridVertex Struct
// Constains information about grid-vertices
#include <iostream>

template<int dim = 0 >
class GridVertex {
public:
    GridVertex();
};

template<> class GridVertex<1> {
    public:
        GridVertex(double pos);
        double position;
        double Efield;
        double Bfield;
        double effectiveCharge;
        
};


#endif 