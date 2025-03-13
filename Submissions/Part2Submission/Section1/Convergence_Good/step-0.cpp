#include <iomanip>

#include "NBodySimulation.h"

/**
 * You can compile this file with
 *   make step-0
 * and run it with
 *   ./step-0
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */
int main (int argc, char** argv) {

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulation nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  nbs.calculateTotalEnergy(true);

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.calculateTotalEnergy(false);

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}

// int main(int argc, char** argv) {
//   NBodySimulation simulation;
//   simulation.setUp(argc, argv);
  
//   // Analyze convergence order
//   // double baseTimeStep = 0.0001;  // Starting time step for analysis
//   // double maxTestTime = 1.0;      // Physical time to simulate for each test
//   // int numRefinements = 5;        // Number of time step refinements to test
//   // double energyThreshold = 0.05; // Allow 5% energy drift before declaring instability
  
//   // simulation.analyzeConvergenceOrder(
//   //     baseTimeStep,
//   //     maxTestTime,
//   //     numRefinements,
//   //     energyThreshold
//   // );
  
//   // // Continue with normal simulation...
//   // simulation.timeStepSize = baseTimeStep; // Use stable time step
  
//   simulation.openParaviewVideoFile();
//   simulation.takeSnapshot();
  
//   while (!simulation.hasReachedEnd()) {
//       simulation.updateBody();
//       simulation.takeSnapshot();
//   }
  
//   simulation.closeParaviewVideoFile();
//   simulation.printSummary();
  
//   return 0;
// }

