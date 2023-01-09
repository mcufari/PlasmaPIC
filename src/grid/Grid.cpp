/*
Implementations of multidimensional grids
See grid.hpp for documentation
*/

#include "Grid.hpp"

Grid<1>::Grid(int NPoints, double dx, double lBoundary) : NP(NPoints), dx(dx), lBound(lBoundary){
    rBound = lBound + (NPoints-1)*dx; //1 point is used at each boundary
    boundaryCondition = "Periodic"; //Only boundary condition permitted at present in periodic

    for(int i = 0; i < NPoints; i++){
        GridVertex<1>* GV = new GridVertex<1>(lBound + i*dx); //Create a new grid vertex object 
        indexVertexMap[i] = GV;  //Insert this object into lookup table
    }

    poissonMatrix = (double*) calloc(NPoints*NPoints, sizeof(double));
    for(int i = 1; i < NPoints-1; i++){
        poissonMatrix[i*NPoints + i] = -2;
        poissonMatrix[i*NPoints + i+1] = 1;
        poissonMatrix[i*NPoints + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    poissonMatrix[NPoints - 1] = 1;

    poissonMatrix[(NPoints*NPoints)-1] = -2;
    poissonMatrix[(NPoints*NPoints)-2] = 1;
    poissonMatrix[(NPoints-1)*NPoints] = 1;
    
    ipiv = (int*) calloc(NPoints, sizeof(int));
    int info = 0;
    int nrhs = 1;
    for(int i = 0; i < NP; i++){
        for(int j = 0; j < NP; j++){
            std::cout << poissonMatrix[i*NP + j] << " ";
        }
        std::cout << std::endl;
    }
    dgetrf(&NP,&NP,poissonMatrix,&NP,ipiv,&info);
    

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
    
    poissonMatrix = (double*) calloc(NPoints*NPoints, sizeof(double));
    for(int i = 1; i < NPoints-1; i++){
        poissonMatrix[i*NPoints + i] = -2;
        poissonMatrix[i*NPoints + i+1] = 1;
        poissonMatrix[i*NPoints + i-1] = 1;
    }
    poissonMatrix[0] = -2;
    poissonMatrix[1] = 1;
    poissonMatrix[NPoints - 1] = 1;

    poissonMatrix[(NPoints*NPoints)-1] = -2;
    poissonMatrix[(NPoints*NPoints)-2] = 1;
    poissonMatrix[(NPoints-1)*NPoints] = 1;
    
    ipiv = (int*) calloc(NPoints, sizeof(int));
    int info = 0;
    int nrhs = 1;
    for(int i = 0; i < NP; i++){
        for(int j = 0; j < NP; j++){
            std::cout << poissonMatrix[i*NP + j] << " ";
        }
        std::cout << std::endl;
    }
    
     dgetrf(&NP,&NP,poissonMatrix,&NP,ipiv,&info);

}

void Grid<1>::vertexInfoTraverse(){
    for(int i = 0; i < NP; i++){
        std::cout << "Pos " << i << " is: "<< indexVertexMap[i]->position << std::endl;
        std::cout << "effectiveCharge " << i << " is: " << indexVertexMap[i]->effectiveCharge << std::endl;
    }
}

void Grid<1>::poissonSolver(){
    const char trans = 'T';
    double* density = (double*) malloc(sizeof(double) * NP);
    for(int i = 0; i < NP; i++){
        density[i] = -1*indexVertexMap[i]->effectiveCharge * dx * dx;
    }
   
    int info = 0;
    int nrhs = 1;
    dgetrf(&NP,&NP,poissonMatrix,&NP,ipiv,&info);
    dgetrs(&trans,&NP,&nrhs,poissonMatrix,&NP,ipiv,density,&NP,&info);

    for(int i = 0; i < NP; i++){
        std::cout << "Phi at " << i << " is: " << density[i] << std::endl;
    }
    free(density);
    
}
