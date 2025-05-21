ROOTDIR := $(shell pwd)
OUTPUTDIR := $(ROOTDIR)/paraview-output/

CXX := nvc++
CXXFLAGS := -mp=gpu -O3 -march=native -std=c++14

TARGETS := step-0 step-1 step-2 step-3 step-4 step-5 step-5-gpu initial_conditions_generator adapted-step-0

.PHONY: all cleanall clean clean_paraview

all: $(TARGETS)

# Create output directory
$(shell mkdir -p $(OUTPUTDIR))

# Object files
OBJS_COMMON := NBodySimulation.o NBodySimulationVectorised.o NBodyInputReader.o
OBJS_STEP5_GPU := NBodySimulationGPU.o step-5-gpu.o NBodyInputReader.o

# Compile rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Build steps
step-0: $(OBJS_COMMON) step-0.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-1: $(OBJS_COMMON) step-1.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-2: $(OBJS_COMMON) step-2.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-3: $(OBJS_COMMON) step-3.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-4: $(OBJS_COMMON) step-4.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-5: step-5.o NBodySimulation.o NBodySimulationInstrumented.o
	$(CXX) $(CXXFLAGS) -o $@ $^

step-5-gpu: $(OBJS_STEP5_GPU)
	$(CXX) $(CXXFLAGS) -o $@ $^

initial_conditions_generator: initial_conditions_generator.o NBodyInputReader.o
	$(CXX) $(CXXFLAGS) -o $@ $^

adapted-step-0: adapted-step-0.o $(OBJS_COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Cleaning rules
cleanall: clean clean_paraview clean_init

clean:
	rm -rf $(ROOTDIR)/step-[0-5] $(ROOTDIR)/step-5-gpu $(ROOTDIR)/initial_conditions_generator \
	       $(ROOTDIR)/adapted-step-0 $(ROOTDIR)/*.o

clean_init:
	rm -rf $(ROOTDIR)/step-[0-5].sh

clean_paraview:
	if test -d "$(OUTPUTDIR)"; then \
		find $(OUTPUTDIR) -iname "result-*.vtp" -delete; \
		rm -f $(OUTPUTDIR)/result.pvd; \
	fi
