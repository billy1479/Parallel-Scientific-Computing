#include <iomanip>
#include "NBodySimulation.h"

/**
 * NBodySimulationCollision implements collision detection and merging
 * with the merge condition |x₁-x₂|≤C(m₁+m₂) where C=10⁻²/N
 */
class NBodySimulationCollision : public NBodySimulation {
private:
  double collisionConstant; // C = 10⁻²/N
  bool* merged; // Array to track which bodies have been merged

public:
  // Constructor initializes collision-specific variables
  NBodySimulationCollision() : NBodySimulation(), collisionConstant(0), merged(nullptr) {}
  
  // Destructor for cleanup
  ~NBodySimulationCollision() {
    if (merged != nullptr) {
      delete[] merged;
    }
  }
  
  // Override setUp to initialize collision parameters
  void setUp(int argc, char** argv) {
    // Call parent class setup first
    NBodySimulation::setUp(argc, argv);
    
    // Initialize the collision constant C = 10^(-2)/N
    collisionConstant = 0.01 / NumberOfBodies;
    
    // Initialize the merged array
    merged = new bool[NumberOfBodies]();  // Initialize all to false
    
    std::cout << "Collision simulation setup complete. Collision constant: " 
              << collisionConstant << std::endl;
  }
  
  // Override updateBody to include collision detection and merging
  void updateBody() {
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    // Calculate forces (same as parent class)
    double* force0 = new double[NumberOfBodies]();
    double* force1 = new double[NumberOfBodies]();
    double* force2 = new double[NumberOfBodies]();

    if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

    for (int i = 0; i < NumberOfBodies; i++) {
      force0[i] = 0.0;
      force1[i] = 0.0;
      force2[i] = 0.0;
    }

    // Calculate forces between pairs of bodies
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i + 1; j < NumberOfBodies; j++) {
        if (i != j) {
          // x,y,z forces acting on particle i
          force0[i] += force_calculation(i, j, 0);
          force1[i] += force_calculation(i, j, 1);
          force2[i] += force_calculation(i, j, 2);
          // x,y,z symmetric forces acting on particle j
          force0[j] -= force_calculation(i, j, 0);
          force1[j] -= force_calculation(i, j, 1);
          force2[j] -= force_calculation(i, j, 2);
        }
      }
    }

    // Update positions first (same as parent class)
    for (int i = 0; i < NumberOfBodies; i++) {
      x[i][0] = x[i][0] + timeStepSize * v[i][0];
      x[i][1] = x[i][1] + timeStepSize * v[i][1];
      x[i][2] = x[i][2] + timeStepSize * v[i][2];
    }

    // Check for collisions and merge bodies if needed
    mergeBodies();

    // Update velocities after merging (adjusted for possibly reduced body count)
    for (int i = 0; i < NumberOfBodies; i++) {
      v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
      v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
      v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

      maxV = std::max(maxV, std::sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]));
    }

    t += timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
  }

private:
  // Helper function to detect and merge colliding bodies
  void mergeBodies() {
    bool mergeOccurred = false;

    for (int i = 0; i < NumberOfBodies; i++) {
      if (merged[i]) continue;  // Skip already merged bodies

      for (int j = i + 1; j < NumberOfBodies; j++) {
        if (merged[j]) continue;  // Skip already merged bodies

        // Calculate the Euclidean distance between bodies i and j
        const double distance = sqrt(
          (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
          (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
          (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])
        );

        // Check if they should merge: |x1-x2| ≤ C(m1+m2)
        if (distance <= collisionConstant * (mass[i] + mass[j])) {
          mergeOccurred = true;

          // Calculate the new position and velocity as the mass-weighted mean
          double totalMass = mass[i] + mass[j];
          
          // Position: x = (m1*x1 + m2*x2)/(m1+m2)
          for (int d = 0; d < 3; d++) {
            x[i][d] = (mass[i] * x[i][d] + mass[j] * x[j][d]) / totalMass;
            
            // Velocity: v = (m1*v1 + m2*v2)/(m1+m2)
            v[i][d] = (mass[i] * v[i][d] + mass[j] * v[j][d]) / totalMass;
          }

          // Update the mass of body i
          mass[i] = totalMass;

          // Mark body j as merged
          merged[j] = true;

          std::cout << "Bodies " << i << " and " << j << " merged with new mass "
                    << mass[i] << " at position ("
                    << x[i][0] << ", " << x[i][1] << ", " << x[i][2] << ")" << std::endl;
        }
      }
    }

    // If merges occurred, we need to compact the arrays
    if (mergeOccurred) {
      int newIndex = 0;
      for (int i = 0; i < NumberOfBodies; i++) {
        if (!merged[i]) {
          // If this body wasn't merged, move it to the next available slot
          if (newIndex != i) {
            for (int d = 0; d < 3; d++) {
              x[newIndex][d] = x[i][d];
              v[newIndex][d] = v[i][d];
            }
            mass[newIndex] = mass[i];
          }
          newIndex++;
        }
      }

      // Update the number of bodies
      int oldNumberOfBodies = NumberOfBodies;
      NumberOfBodies = newIndex;

      std::cout << "After merging, number of bodies reduced from "
                << oldNumberOfBodies << " to " << NumberOfBodies << std::endl;

      // Reset the merged array for future use
      for (int i = 0; i < oldNumberOfBodies; i++) {
        merged[i] = false;
      }
    }
  }
};

/**
 * Main routine.
 */
int main(int argc, char** argv) {
  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationCollision nbs;
  nbs.setUp(argc, argv);
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