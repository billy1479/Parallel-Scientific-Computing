/**
 * You can compile this file with
 *   make step-4
 * and run it with
 *   ./step-4
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

 #include <iomanip>
 #include <vector>
 #include <memory>
 #include <cmath>
 #include <iostream>
 #include <limits>
 #include "NBodySimulation.h"
 #include "OctreeNode.h"
 
 /**
  * Implementation of the NBodySimulationCollision class from step-3
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
     collisionConstant = 0.0001 / NumberOfBodies;
     
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
 
 protected:
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
  * Barnes-Hut simulation for N-body problems
  * Inherits collision detection functionality from NBodySimulationCollision
  */
 class NBodyBarnsHutSimulation : public NBodySimulationCollision {
 private:
     OctreeNode* octree;       // The octree for Barnes-Hut algorithm
     double domainSize;        // Size of the simulation domain
     
 public:
     // Constructor
     NBodyBarnsHutSimulation() : NBodySimulationCollision(), octree(nullptr), domainSize(0.0) {}
     
     // Destructor
     ~NBodyBarnsHutSimulation() {
         if (octree != nullptr) {
             delete octree;
             octree = nullptr;
         }
     }
     
     // Override setUp to initialize Barnes-Hut specific parameters
     void setUp(int argc, char** argv) {
         // Call parent class setup first
         NBodySimulationCollision::setUp(argc, argv);
         
         // Initialize Barnes-Hut specific members
         octree = nullptr;
         domainSize = 0.0;
         
         std::cout << "Barnes-Hut simulation initialized" << std::endl;
     }
     
     // Override updateBody to use Barnes-Hut algorithm
     void updateBody() {
         try {
             // Initial setup
             timeStepCounter++;
             maxV = 0.0;
             minDx = std::numeric_limits<double>::max();
 
             // Build the octree for this time step
             buildOctree();
 
             // Allocate and initialize force arrays
             double* force0 = new double[NumberOfBodies]();
             double* force1 = new double[NumberOfBodies]();
             double* force2 = new double[NumberOfBodies]();
 
             // Reset forces
             for (int i = 0; i < NumberOfBodies; i++) {
                 force0[i] = 0.0;
                 force1[i] = 0.0;
                 force2[i] = 0.0;
             }
 
             // Compute forces for each body using Barnes-Hut algorithm
             for (int i = 0; i < NumberOfBodies; i++) {
                 // Skip bodies with invalid positions
                 bool invalidPosition = false;
                 for (int d = 0; d < 3; d++) {
                     if (std::isnan(x[i][d]) || std::isinf(x[i][d])) {
                         invalidPosition = true;
                         break;
                     }
                 }
                 if (invalidPosition) continue;
                 
                 // Use a temporary array for this body's forces
                 double forces[3] = {0.0, 0.0, 0.0};
                 
                 // Compute force using the octree
                 if (octree != nullptr) {
                     octree->computeForce(i, x, forces, mass);
                 }
                 
                 // Store the forces
                 force0[i] = forces[0];
                 force1[i] = forces[1];
                 force2[i] = forces[2];
             }
 
             // Update positions using current velocities
             for (int i = 0; i < NumberOfBodies; i++) {
                 x[i][0] = x[i][0] + timeStepSize * v[i][0];
                 x[i][1] = x[i][1] + timeStepSize * v[i][1];
                 x[i][2] = x[i][2] + timeStepSize * v[i][2];
             }
 
             // Check for collisions and merge bodies if needed
             mergeBodies();
 
             // Update velocities using computed forces
             for (int i = 0; i < NumberOfBodies; i++) {
                 v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
                 v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
                 v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];
 
                 // Update maximum velocity
                 maxV = std::max(maxV, std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]));
             }
 
             // Recalculate minimum distance between bodies (for reporting purposes)
             if (NumberOfBodies > 1) {
                 for (int i = 0; i < NumberOfBodies - 1; i++) {
                     for (int j = i + 1; j < NumberOfBodies; j++) {
                         double distance = std::sqrt(
                             (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                             (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                             (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                         );
                         minDx = std::min(minDx, distance);
                     }
                 }
             } else {
                 minDx = 0.0;  // No distances to calculate for a single body
             }
 
             // Update time
             t += timeStepSize;
             
             // Clean up
             delete[] force0;
             delete[] force1;
             delete[] force2;
             
         } catch (std::exception& e) {
             std::cerr << "Exception in updateBody: " << e.what() << std::endl;
         } catch (...) {
             std::cerr << "Unknown exception in updateBody" << std::endl;
         }
     }
 
 private:
     // Build the Barnes-Hut octree
     void buildOctree() {
         try {
             // Delete old octree if it exists
             if (octree != nullptr) {
                 delete octree;
                 octree = nullptr;
             }
             
             // Calculate the domain size
             domainSize = findDomainSize();
             
             // Create new octree
             octree = new OctreeNode(0.0, 0.0, 0.0, domainSize / 2.0);
             
             // Insert all bodies into the octree
             for (int i = 0; i < NumberOfBodies; i++) {
                 // Check for invalid positions
                 bool invalidPosition = false;
                 for (int d = 0; d < 3; d++) {
                     if (std::isnan(x[i][d]) || std::isinf(x[i][d])) {
                         invalidPosition = true;
                         break;
                     }
                 }
                 
                 // Skip bodies with invalid positions
                 if (invalidPosition) continue;
                 
                 // Insert into the octree
                 octree->insert(i, x, mass);
             }
             
         } catch (std::exception& e) {
             std::cerr << "Exception in buildOctree: " << e.what() << std::endl;
         } catch (...) {
             std::cerr << "Unknown exception in buildOctree" << std::endl;
         }
     }
     
     // Calculate appropriate domain size for the octree
     double findDomainSize() {
         try {
             // Start with a reasonable default size
             double maxCoord = 1.0;
             
             // Find maximum absolute coordinate among all particles
             for (int i = 0; i < NumberOfBodies; i++) {
                 for (int d = 0; d < 3; d++) {
                     if (std::isnan(x[i][d]) || std::isinf(x[i][d])) {
                         continue;  // Skip invalid coordinates
                     }
                     maxCoord = std::max(maxCoord, std::abs(x[i][d]));
                 }
             }
             
             // Add margin and ensure it's reasonable
             double size = maxCoord * 2.0 * 1.2;  // 2x for diameter, 20% margin
             
             // Limit to reasonable size
             size = std::max(1.0, std::min(size, 1000.0));
             
             // Round up to power of 2
             double power = std::ceil(std::log2(size));
             return std::pow(2.0, power);
         } catch (std::exception& e) {
             std::cerr << "Exception in findDomainSize: " << e.what() << std::endl;
             return 100.0;  // Reasonable default
         }
     }
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
  NBodyBarnsHutSimulation nbs;
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
