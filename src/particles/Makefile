CXX = icc
CXXFLAGS = -g -qmkl -fopenmp -O0 -march=core-avx2

Particle.o : Particle.cpp Particle.hpp
	$(CXX) -o $@ Particle.cpp $(CXXFLAGS) 

particleInCell.o: particleInCell.cpp particleInCell.hpp
	$(CXX) -o $@ particleInCell.cpp $(CXXFLAGS)

