# N-body Simulation Project

This repository contains the implementation of an N-body simulation system for the Parallel Scientific Computing (COMP3741) course. The project is divided into three main components, each in its own directory:

## Project Structure

```
.
├── README.md
├── parallel_computing/     # First component - due Jan 21, 2025
├── numerical_algorithms/   # Second component - due Mar 13, 2025
└── system_administration/  # Third component - due Apr 29, 2025
```

## 1. Parallel Computing

This folder contains the implementation of vectorized and parallelized versions of the N-body simulation:

- `step-0.cpp` - Original non-vectorized, non-parallelized baseline code
- `step-1.cpp` - Vectorized implementation using OpenMP SIMD directives
- `step-2.cpp` - Parallelized implementation using OpenMP multi-threading
- `NBodySimulation.h/cpp` - Base simulation class
- `NBodySimulationVectorised.cpp` - Vectorized simulation class
- `NBodySimulationParallelised.cpp` - Parallelized simulation class
- `create_initial_conditions.py` - Script for generating initial configurations
- `Makefile` - Compilation instructions

### Usage

```bash
# Compile all versions
make all

# Run baseline simulation
./step-0 <number_of_bodies> <number_of_timesteps>

# Run vectorized simulation
./step-1 <number_of_bodies> <number_of_timesteps>

# Run parallelized simulation
./step-2 <number_of_bodies> <number_of_timesteps>
```

## 2. Numerical Algorithms

This folder contains implementations focusing on the theoretical aspects and extensions of the N-body solver:

- `step-3.cpp` - Improved N-body solver with collision detection and handling
- `step-4.cpp` - Advanced N-body solver with additional optimizations
- `NBodySimulationCollision.cpp` - Simulation class with collision handling
- `NBodySimulationAdvanced.cpp` - Simulation class with advanced methods
- `Makefile` - Compilation instructions

### Usage

```bash
# Compile all versions
make all

# Run simulation with collision detection
./step-3 <number_of_bodies> <number_of_timesteps>

# Run advanced simulation
./step-4 <number_of_bodies> <number_of_timesteps>
```

## 3. System Administration

This folder contains GPU implementation and cluster computing scripts:

- `step-5-gpu.cpp` - GPU-accelerated N-body simulation using OpenMP
- `install_numactl.sh` - Script to build and install numactl from source
- `run_all.sh` - Batch script to generate and run multiple test cases
- `Makefile` - Compilation instructions with GPU targets

### Usage

```bash
# Compile CPU version
make step-0

# Compile GPU version
make step-5-gpu

# Install numactl locally
./install_numactl.sh

# Run batch tests
sbatch run_all.sh
```

## Running on Durham's Supercomputers

All code is designed to run on Hamilton (for parallel computing and numerical algorithms components) and NCC (for system administration component).

### Hamilton

Example job submission:

```bash
# Submit a job
sbatch --nodes=1 --ntasks-per-node=16 --time=00:30:00 --partition=test ./run_job.sh
```

### NCC

GPU execution example:

```bash
# Run on GPU
./step-5-gpu <number_of_bodies> <number_of_timesteps>

# Verify numactl architecture information
numactl -H
```

## Notes

- Make sure to follow the assignment's submission requirements when submitting your work.
- Each component requires a separate report (max 1 page excluding figures) answering specific questions.
- The code must compile and run on Hamilton/NCC using the provided Makefile.