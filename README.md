# Parallel-Computing-Project
This work incorporates an MPI parallel Laplace equation solver for high-performance computing environments with OpenMP and serial comparison. The solver was built to examine performance trends of various parallelization methods on 2D finite difference problems.

This project implements and compares three different approaches to solving the 2D Laplace equation:

- Serial Implementation: Single-threaded baseline version
- OpenMP Implementation: Shared-memory parallelization using OpenMP
- MPI Implementation: Distributed-memory parallelization using Message Passing Interface

The solver uses the Jacobi iterative method with a 5-point stencil to solve the steady-state heat distribution problem on a 2D grid.

Key Findings:

- OpenMP demonstrates superior performance over MPI across all matrix sizes
- MPI shows poor scalability due to communication overhead
- OpenMP achieves best speedup (1.43x) on smaller matrices
- Performance gap between OpenMP and MPI increases with matrix size

Required Software

- GCC Compiler: Version 4.4.7 (The supercomputer used from SINES lab had this compiler so my code is optimised for this compiler)
- OpenMP: Usually included with GCC

# System Requirements

Linux/Unix environment
Multi-core CPU for parallel execution

# Installation & Build:

# Clone the repository:
git clone https://github.com/princess-humario/Parallel-Computing-Project
cd Parallel-Computing-Project

# Build all versions:
make clean
make all

# Usage:
# Serial Version
./LaplaceSolver
//Input: rows, columns, max_iterations, tolerance
# OpenMP Version
./LaplaceSolver_super
//Input: rows, columns, max_iterations, tolerance, num_threads
# MPI Version
mpirun -np <num_processes> ./LaplaceSolver_mpi
//Input: rows, columns, max_iterations, tolerance
