#ifndef NBODY_SIMULATION_H
#define NBODY_SIMULATION_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

// Forward declaration of SimulationParams
struct SimulationParams;

class NBodySimulation {

  public:
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
  double*  mass;

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

  NBodySimulation ();
  ~NBodySimulation ();

  void checkInput(int argc, char** argv);
  void setUp (int argc, char** argv);
  void setUpFromParams(const SimulationParams& params);
  double force_calculation (int i, int j, int direction);
  void updateBody ();
  bool hasReachedEnd ();
  void takeSnapshot ();
  void openParaviewVideoFile ();
  void closeParaviewVideoFile ();
  void printParaviewSnapshot ();
  void printSnapshotSummary ();
  void printSummary ();

};

#endif // NBODY_SIMULATION_H