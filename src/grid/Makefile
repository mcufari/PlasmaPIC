CXX = icc
CXXFLAGS = -g -qmkl -fopenmp -O3 -march=core-avx2

Grid.o : Grid.cpp Grid.hpp
	$(CXX) -o $@ Grid.cpp $(CXXFLAGS)

