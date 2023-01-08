CXX = icc
CXXFLAGS = -g -qmkl

Particle.o : Particle.cpp Particle.hpp
	$(CXX) -o $@ Particle.cpp $(CXXFLAGS) 

particleInCell.o: particleInCell.cpp particleInCell.hpp
	$(CXX) -o $@ particleInCell.cpp $(CXXFLAGS)

