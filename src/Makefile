#src

CXX = icc
CXXFLAGS = -g -qmkl -fopenmp -O3 -march=core-avx2
OBJECTS = main.o

main.o: main.cpp
	$(CXX) -o $@ main.cpp $(CXXFLAGS) 
