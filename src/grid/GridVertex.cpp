#include "GridVertex.hpp"


GridVertex<1>::GridVertex(double pos): position(pos) {
    this->Bfield=0;
    this->Efield=0;
    this->effectiveCharge=0;
}



