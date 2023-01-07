#include "Grid.hpp"

Grid<1>::Grid(int NPoints, double dx, double lBoundary) : NP(NPoints), dx(dx), lBound(lBoundary){
    rBound = lBound + (NPoints-1)*dx; //1 point is used at each boundary
    boundaryCondition = "Periodic";

    for(int i = 0; i < NPoints; i++){
        GridVertex<1>* GV = new GridVertex<1>(lBound + i*dx);
        indexVertexMap[i] = GV; 
    }
}



Grid<1>::Grid(int NPoints, double dx, double lBoundary, std::string boundCondition) : 
NP(NPoints), 
dx(dx), 
lBound(lBoundary), 
boundaryCondition(boundCondition) 
{
    rBound = lBound + (NPoints - 1) * dx;
    for(int i = 0; i < NPoints; i++){
        GridVertex<1>* GV = new GridVertex<1>(lBound + i*dx);
        indexVertexMap[i] = GV; 
    }
}

void Grid<1>::vertexInfoTraverse(){
    for(int i = 0; i < NP; i++){
        std::cout << "Pos " << i << " is: "<< indexVertexMap[i]->position << std::endl;
        std::cout << "effectiveCharge " << i << " is: " << indexVertexMap[i]->effectiveCharge << std::endl;
    }
}
