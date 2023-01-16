CXX = icc

OBJECTS = src/particles/Particle.o src/particles/particleInCell.o src/grid/Grid.o src/main.o

CXXFLAGS = -g -qmkl -fopenmp -O3 -march=core-avx2

plasma: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CXXFLAGS)

clean:
	rm -f $(OBJECTS) plasma