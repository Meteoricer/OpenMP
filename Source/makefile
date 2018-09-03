CC=g++
CFLAGS=-fopenmp -O0 -std=c++11 -c -Wall
LDFLAGS=-fopenmp -O0
SOURCES=main.cpp Anchor_Unit.cpp DeconstructionFunction.cpp Diffusion_Unit.cpp General_Functions.cpp Grid_Unit.cpp initialize.cpp Polymer_Bundle.cpp Polymer_Cluster.cpp Polymer_Unit.cpp Rule_Anchoring_Deanchoring.cpp Rule_Bundling_Debundling.cpp Rule_Diffusion.cpp Rule_Fragmentation_Annealing.cpp Rule_Hydrolysis.cpp Rule_Polymerize_Depolymerize.cpp Rule_Rotation.cpp Rule_Structure.cpp SaveandLoad.cpp test.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=test

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
