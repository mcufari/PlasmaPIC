#ifndef GRID_HEADER
#define GRID_HEADER
#include "GridVertex.hpp"
#include <iostream>
#include <string>
#include <unordered_map>

template<int dim = 0>
class Grid{
public:
    Grid(); 
};


template<> class Grid<1>{
    public:
    Grid(int NPoints, double dx, double lBoundary);
 
    Grid(int NPoints, double dx, double lBoundary, std::string boundCondition);
  
    void vertexInfoTraverse();
    
    std::unordered_map<int, GridVertex<1>*> indexVertexMap;

    int NP;
    double dx;
    double lBound;
    double rBound;
    std::string boundaryCondition;

    
};


#endif 