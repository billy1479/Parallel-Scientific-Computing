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

//   while (!nbs.hasReachedEnd()) {
//     nbs.updateBody();
//     nbs.takeSnapshot();
//   }

//   nbs.printSummary();
//   nbs.closeParaviewVideoFile();

//   return 0;
// }


#include <vector>
#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    // Define different time steps for testing
    std::vector<double> dt_values = {0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 0.0009765625, 0.00048828125, 0.000244140625};  // Different time steps
    std::vector<double> errors;
    
    // Run a reference simulation with a very small dt
    double dt_ref = 1;
    NBodySimulation refSim;
    refSim.setUp(argc, argv, dt_ref);  // Setup reference simulation
    refSim.openParaviewVideoFile();
    refSim.takeSnapshot();
    while (!refSim.hasReachedEnd()) {
      refSim.updateBody();
      refSim.takeSnapshot();
    }

    // Save the reference positions
    double** reference_x = refSim.x;

    // Now run simulations for each of the dt values and compute the error
    for (double dt : dt_values) {
        NBodySimulation sim;
        sim.setUp(argc, argv, dt);  // Setup simulation with current dt
        sim.openParaviewVideoFile();
        sim.takeSnapshot();
        while (!sim.hasReachedEnd()) {
          sim.updateBody();
          sim.takeSnapshot();
        } // Run the simulation to final time
        
        // Compute error compared to the reference solution
        double error = sim.computeError(sim.x, reference_x, sim.NumberOfBodies);
        errors.push_back(error);
        
        std::cout << "Error for dt=" << dt << " is " << error << std::endl;
    }

    // Compute and display convergence order
    for (size_t i = 1; i < dt_values.size(); i++) {
        double order = log(errors[i] / errors[i - 1]) / log(dt_values[i] / dt_values[i - 1]);
        std::cout << "Convergence order between dt=" << dt_values[i-1]
                  << " and dt=" << dt_values[i] << " is " << order << std::endl;
    }

    return 0;
}
