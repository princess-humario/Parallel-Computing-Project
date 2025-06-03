CXX = g++
MPICXX = mpic++
CXXFLAGS = -Wall -O2
OPENMP_FLAGS = -fopenmp

# Targets
all: serial openmp mpi

serial: src/LaplaceSolver.cpp
	$(CXX) $(CXXFLAGS) -o LaplaceSolver src/LaplaceSolver.cpp

openmp: src/LaplaceSolver_super.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) -o LaplaceSolver_super src/LaplaceSolver_super.cpp

mpi: src/mpi.cpp
	$(MPICXX) $(CXXFLAGS) -o LaplaceSolver_mpi src/mpi.cpp

clean:
	rm -f LaplaceSolver LaplaceSolver_super LaplaceSolver_mpi
