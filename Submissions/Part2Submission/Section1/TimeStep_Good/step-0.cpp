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
// int main (int argc, char** argv) {

//   std::cout << std::setprecision(15);

//   // Code that initialises and runs the simulation.
//   NBodySimulation nbs;
//   nbs.setUp(argc,argv);
//   nbs.openParaviewVideoFile();
//   nbs.takeSnapshot();

//   nbs.calculateTotalEnergy(true);

//   while (!nbs.hasReachedEnd()) {
//     nbs.updateBody();
//     nbs.takeSnapshot();
//   }

//   nbs.calculateTotalEnergy(false);

//   nbs.printSummary();
//   nbs.closeParaviewVideoFile();

//   return 0;
// }

int main(int argc, char** argv) {
  NBodySimulation simulation;
  simulation.setUp(argc, argv);
  
  // Find largest stable time step
  // double initialTimeStep = 0.00001;  // Start with a very small time step
  // double maxTimeStep = 1;          // Maximum time step to test
  // double maxTestTime = 60.0;          // Physical time to simulate for stability test
  // double energyThreshold = 0.01;     // Allow 1% energy drift
  
  // double stableTimeStep = simulation.findLargestStableTimeStep(
  //     initialTimeStep, 
  //     maxTimeStep, 
  //     maxTestTime, 
  //     energyThreshold
  // );
  
  // // Use the found stable time step
  // simulation.timeStepSize = stableTimeStep;
  
  // Now run with the determined step size
  simulation.openParaviewVideoFile();
  simulation.takeSnapshot();

  simulation.calculateTotalEnergy(true);
  
  while (!simulation.hasReachedEnd()) {
      simulation.updateBody();
      simulation.takeSnapshot();
  }

  simulation.calculateTotalEnergy(false);
  
  simulation.closeParaviewVideoFile();
  simulation.printSummary();
  
  return 0;
}

