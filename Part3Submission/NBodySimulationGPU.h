#ifndef NBODY_SIMULATION_GPU_H
#define NBODY_SIMULATION_GPU_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

class NBodySimulationGPU {

 private:
  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies;

  /**
   * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
   * each pointer represents one molecule/particle/body.
   */
  double** x;

  /**
   * Equivalent to x storing the velocities.
   */
  double** v;

  /**
   * One mass entry per molecule/particle.
   */
  double* mass;

  /**
   * Force arrays for CPU
   */
  double* force0;
  double* force1;
  double* force2;

  /**
   * Flat arrays for GPU computation:
   * x_flat: [x1,y1,z1,x2,y2,z2,...] - positions
   * v_flat: [vx1,vy1,vz1,vx2,vy2,vz2,...] - velocities
   * force_flat: [fx1,fy1,fz1,fx2,fy2,fz2,...] - forces
   */
  double* x_flat;
  double* v_flat;
  double* force_flat;

  /**
   * Global time step size used.
   */
  double timeStepSize;

  /**
   * Maximum velocity of all particles.
   */
  double maxV;

  /**
   * Minimum distance between two elements.
   */
  double minDx;

  /**
   * Stream for video output file.
   */
  std::ofstream videoFile;

  /**
   * Output counters.
   */
  int snapshotCounter;
  int timeStepCounter;

  /**
   * Allocate memory for force arrays
   */
  void allocateForceArrays();

  /**
   * Free force arrays
   */
  void freeForceArrays();

  /**
   * Allocate and initialize GPU flat arrays
   */
  void setupGPUArrays();

  /**
   * Free GPU arrays
   */
  void freeGPUArrays();

  /**
   * Convert from CPU's 2D arrays to GPU's flat arrays
   */
  void convertToFlatArrays();

  /**
   * Convert from GPU's flat arrays back to CPU's 2D arrays
   */
  void convertFromFlatArrays();

 public:
  NBodySimulationGPU();
  ~NBodySimulationGPU();

  /**
   * Check that the number command line parameters is correct.
   */
  void checkInput(int argc, char** argv);

  /**
   * Set up scenario from the command line or input file.
   */
  void setUp(int argc, char** argv);
  void setupFromCommandLine(int argc, char** argv);
  void setupFromInputFile(const std::string& filename);

  /**
   * Compute forces (CPU version - kept for reference).
   */
  double force_calculation(int i, int j, int direction);

  /**
   * Run the simulation step using GPU acceleration.
   */
  void updateBody();
  void runSimulation();

  /**
   * Check if the last time step has been reached (simulation is completed).
   */
  bool hasReachedEnd();

  /**
   * Take simulation snapshots and print summary to standard output.
   */
  void takeSnapshot();

  /**
   * Handle Paraview output.
   */
  void openParaviewVideoFile();
  void closeParaviewVideoFile();
  void printParaviewSnapshot();

  /**
   * Handle terminal output.
   */
  void printSnapshotSummary();
  void printSummary();
};

#endif // NBODY_SIMULATION_GPU_H