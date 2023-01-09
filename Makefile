CXX = icc

OBJECTS = src/particles/Particle.o src/particles/particleInCell.o src/grid/Grid.o src/main.o

CXXFLAGS = -g -qmkl

plasma: $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(CXXFLAGS)

clean:
	rm -f $(OBJECTS) plasma