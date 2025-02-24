#include <iomanip>
#include "NBodySimulation.h"

/**
 * You can compile this file with
 *   make step-3.
 * and run it with
 *   ./step-3
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * Class that extends NBodySimulation to handle object collisions.
 * Objects merge when they collide, with mass conservation and
 * mass-weighted averaging of positions and velocities.
 */
class NBodySimulationCollision : public NBodySimulation {
private:
    int initialNumberOfBodies; // Store the initial number of bodies for collision tolerance
    double collisionTolerance; // The collision distance factor C

    /**
     * Checks for collisions between all pairs of bodies and merges them if necessary.
     * @return Number of collisions detected and processed
     */
    int handleCollisions() {
        if (NumberOfBodies <= 1) return 0; // No collisions possible with 0 or 1 bodies
        
        int collisionCount = 0;
        bool collisionOccurred = true;
        
        // Keep checking for collisions until no more are found in a pass
        // This handles the case where after merging A and B, the resulting body
        // might immediately collide with body C
        while (collisionOccurred) {
            collisionOccurred = false;
            
            // Create a list of collision pairs to process
            std::vector<std::pair<int, int>> collisionPairs;
            
            // Check all pairs of bodies for collisions
            for (int i = 0; i < NumberOfBodies; i++) {
                for (int j = i + 1; j < NumberOfBodies; j++) {
                    // Calculate distance between bodies
                    double dx = x[i][0] - x[j][0];
                    double dy = x[i][1] - x[j][1];
                    double dz = x[i][2] - x[j][2];
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    
                    // Calculate collision threshold based on masses
                    double collisionThreshold = collisionTolerance * (mass[i] + mass[j]);
                    
                    // Check if collision occurred
                    if (distance <= collisionThreshold) {
                        collisionPairs.push_back(std::make_pair(i, j));
                        collisionOccurred = true;
                    }
                }
            }
            
            // Process each collision
            for (const auto& pair : collisionPairs) {
                int i = pair.first;
                int j = pair.second;
                
                // Only process if both indices are still valid
                // (they might have been removed in a previous collision)
                if (i < NumberOfBodies && j < NumberOfBodies) {
                    mergeObjects(i, j);
                    collisionCount++;
                }
            }
        }
        
        return collisionCount;
    }

    /**
     * Merges two objects based on conservation of mass and momentum.
     * Object j is merged into object i, and object j is removed.
     * 
     * @param i Index of the first object
     * @param j Index of the second object to merge into the first
     */
    void mergeObjects(int i, int j) {
        // Safety check
        if (i >= NumberOfBodies || j >= NumberOfBodies || i == j) {
            return;
        }
        
        // Calculate total mass
        double totalMass = mass[i] + mass[j];
        
        // Calculate mass-weighted position
        for (int d = 0; d < 3; d++) {
            x[i][d] = (mass[i] * x[i][d] + mass[j] * x[j][d]) / totalMass;
        }
        
        // Calculate mass-weighted velocity
        for (int d = 0; d < 3; d++) {
            v[i][d] = (mass[i] * v[i][d] + mass[j] * v[j][d]) / totalMass;
        }
        
        // Update mass of merged object
        mass[i] = totalMass;
        
        // Remove the second object by shifting all objects after it
        removeBody(j);
        
        // Note: NumberOfBodies is reduced by 1 in removeBody()
    }
    
    /**
     * Removes a body from the simulation by shifting all subsequent bodies.
     * 
     * @param index Index of the body to remove
     */
    void removeBody(int index) {
        if (index >= NumberOfBodies) return;
        
        // Shift positions
        for (int i = index; i < NumberOfBodies - 1; i++) {
            for (int d = 0; d < 3; d++) {
                x[i][d] = x[i+1][d];
                v[i][d] = v[i+1][d];
            }
            mass[i] = mass[i+1];
        }
        
        // Decrement body count
        NumberOfBodies--;
    }

public:
    /**
     * Constructor.
     */
    NBodySimulationCollision() : NBodySimulation() {
        initialNumberOfBodies = 0;
        collisionTolerance = 0.0;
    }
    
    /**
     * Sets up the simulation with command line arguments.
     * This overrides the base class method to initialize collision parameters.
     */
    void setUp(int argc, char** argv) {
        // Call the parent class's setUp method
        NBodySimulation::setUp(argc, argv);
        
        // Store the initial number of bodies for collision tolerance calculation
        initialNumberOfBodies = NumberOfBodies;
        
        // Calculate collision tolerance C = 10^(-2) / N
        collisionTolerance = 1e-2 / initialNumberOfBodies;
        
        std::cout << "Collision simulation initialized with:" << std::endl;
        std::cout << "  Initial number of bodies: " << initialNumberOfBodies << std::endl;
        std::cout << "  Collision tolerance C: " << collisionTolerance << std::endl;
    }
    
    /**
     * Updates body positions and velocities for a single time step,
     * and handles collisions after the update.
     * This overrides the base class method to add collision handling.
     */
    void updateBody() {
        // First, update positions and velocities using the base class method
        NBodySimulation::updateBody();
        
        // Then, handle collisions
        int collisions = handleCollisions();
        
        // If any collisions occurred, log them
        if (collisions > 0) {
            std::cout << "Step " << timeStepCounter 
                      << ": " << collisions << " collision(s) occurred. "
                      << "Bodies remaining: " << NumberOfBodies << std::endl;
        }
    }
    
    /**
     * Prints additional collision-related summary information.
     */
    void printSummary() {
        NBodySimulation::printSummary();
        
        std::cout << "Collision Summary:" << std::endl;
        std::cout << "Initial number of bodies: " << initialNumberOfBodies << std::endl;
        std::cout << "Final number of bodies: " << NumberOfBodies << std::endl;
        std::cout << "Number of bodies merged: " << (initialNumberOfBodies - NumberOfBodies) << std::endl;
        std::cout << "Collision tolerance used: " << collisionTolerance << std::endl;
    }
};

class NBodySimulationMolecularForces : public NBodySimulation {

};

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
  NBodySimulationCollision nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
