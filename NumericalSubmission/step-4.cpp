#include <iostream>
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>
#include <immintrin.h>
#include "OctreeNode.h"

/**
 * Memory-aligned structure for better cache performance
 */
struct alignas(64) Body {
  double x, y, z;        // Position
  double vx, vy, vz;     // Velocity
  double mass;           // Mass
  double padding;        // Maintain 64-byte alignment
};

/**
 * NBodySimulationCollision implements collision detection and merging
 * with the merge condition |x₁-x₂|≤C(m₁+m₂) where C=10⁻²/N
 */
class NBodySimulationCollision {
protected:
  Body* bodies;           // Memory-aligned array of bodies
  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;
  int NumberOfBodies;
  int timeStepCounter;
  double timeStepSize;
  double maxV;
  double minDx;
  std::ofstream videoFile;
  int snapshotCounter;
  
  double collisionConstant; // C = 10⁻²/N
  bool* merged; // Array to track which bodies have been merged

public:
  // Constructor
  NBodySimulationCollision() :
      t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
      timeStepCounter(0), timeStepSize(0), maxV(0), minDx(0),
      snapshotCounter(0), collisionConstant(0), merged(nullptr) {
          bodies = nullptr;
      }
  
  // Destructor
  virtual ~NBodySimulationCollision() {
    if (bodies) {
        _mm_free(bodies);  // Free aligned memory
    }
    if (merged != nullptr) {
      delete[] merged;
    }
  }
  
  void setUp(int argc, char** argv) {
    checkInput(argc, argv);

    NumberOfBodies = (argc-4) / 7;

    // Allocate aligned memory for bodies
    bodies = static_cast<Body*>(_mm_malloc(NumberOfBodies * sizeof(Body), 64));
    if (!bodies) {
        std::cerr << "Failed to allocate aligned memory" << std::endl;
        exit(-1);
    }

    int readArgument = 1;

    tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
    tFinal       = std::stof(argv[readArgument]); readArgument++;
    timeStepSize = std::stof(argv[readArgument]); readArgument++;

    for (int i=0; i<NumberOfBodies; i++) {
        bodies[i].x = std::stof(argv[readArgument]); readArgument++;
        bodies[i].y = std::stof(argv[readArgument]); readArgument++;
        bodies[i].z = std::stof(argv[readArgument]); readArgument++;

        bodies[i].vx = std::stof(argv[readArgument]); readArgument++;
        bodies[i].vy = std::stof(argv[readArgument]); readArgument++;
        bodies[i].vz = std::stof(argv[readArgument]); readArgument++;

        bodies[i].mass = std::stof(argv[readArgument]); readArgument++;

        if (bodies[i].mass <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta<=0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    }
    else {
        std::cout << "plot initial setup plus every " << tPlotDelta
                << " time units" << std::endl;
        tPlot = 0.0;
    }
    
    // Initialize collision parameters
    collisionConstant = 0.01 / NumberOfBodies;
    merged = new bool[NumberOfBodies]();  // Initialize all to false
    
    std::cout << "Collision simulation setup complete. Collision constant: " 
              << collisionConstant << std::endl;
  }
  
  // Update method with direct force calculation and collision detection
  virtual void updateBody() {
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    // Force arrays
    double* force0 = new double[NumberOfBodies]();
    double* force1 = new double[NumberOfBodies]();
    double* force2 = new double[NumberOfBodies]();

    if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

    // Calculate forces between pairs of bodies
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i + 1; j < NumberOfBodies; j++) {
        // Calculate distance components once for all three dimensions
        double dx = bodies[j].x - bodies[i].x;
        double dy = bodies[j].y - bodies[i].y;
        double dz = bodies[j].z - bodies[i].z;
        
        // Calculate squared distance
        double distSqr = dx*dx + dy*dy + dz*dz;
        double dist = std::sqrt(distSqr);
        
        // Update minimum distance
        minDx = std::min(minDx, dist);
        
        // Calculate force magnitude once
        double forceMagnitude = bodies[i].mass * bodies[j].mass / (dist * distSqr);
        
        // Calculate force components for all three dimensions at once
        double fx = dx * forceMagnitude;
        double fy = dy * forceMagnitude;
        double fz = dz * forceMagnitude;
        
        // Apply forces to both bodies symmetrically (Newton's third law)
        force0[i] += fx;
        force1[i] += fy;
        force2[i] += fz;
        
        // Apply equal and opposite forces to the other body
        force0[j] -= fx;
        force1[j] -= fy;
        force2[j] -= fz;
      }
    }

    // Update positions first
    for (int i = 0; i < NumberOfBodies; i++) {
      bodies[i].x += timeStepSize * bodies[i].vx;
      bodies[i].y += timeStepSize * bodies[i].vy;
      bodies[i].z += timeStepSize * bodies[i].vz;
    }

    // Check for collisions and merge bodies if needed
    mergeBodies();

    // Update velocities after merging
    for (int i = 0; i < NumberOfBodies; i++) {
      bodies[i].vx += timeStepSize * force0[i] / bodies[i].mass;
      bodies[i].vy += timeStepSize * force1[i] / bodies[i].mass;
      bodies[i].vz += timeStepSize * force2[i] / bodies[i].mass;

      maxV = std::max(maxV, std::sqrt(bodies[i].vx * bodies[i].vx + 
                                      bodies[i].vy * bodies[i].vy + 
                                      bodies[i].vz * bodies[i].vz));
    }

    t += timeStepSize;

    delete[] force0;
    delete[] force1;
    delete[] force2;
  }

  bool hasReachedEnd() {
    return t > tFinal;
  }

  void takeSnapshot() {
    if (t >= tPlot) {
      printParaviewSnapshot();
      printSnapshotSummary();
      tPlot += tPlotDelta;
    }
  }

  void openParaviewVideoFile() {
    videoFile.open("paraview-output/result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
            << " version=\"0.1\""
            << " byte_order=\"LittleEndian\""
            << " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
  }

  void closeParaviewVideoFile() {
    videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
    videoFile.close();
  }

  void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename, filename_nofolder;
    filename << "paraview-output/result-" << counter << ".vtp";
    filename_nofolder << "result-" << counter << ".vtp";
    std::ofstream out(filename.str().c_str());
    
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float64\""
        << " NumberOfComponents=\"3\""
        << " format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++) {
        out << bodies[i].x << " "
            << bodies[i].y << " "
            << bodies[i].z << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    out.close();

    videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
  }

  void printSnapshotSummary() {
    std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t=" << t
            << ",\t dt=" << timeStepSize
            << ",\t v_max=" << maxV
            << ",\t dx_min=" << minDx
            << std::endl;
  }

  void printSummary() {
    std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
    std::cout << "Position of first remaining object: "
            << bodies[0].x << ", " << bodies[0].y << ", " << bodies[0].z << std::endl;
  }

  void checkInput(int argc, char** argv) {
    if (argc==1) {
      std::cerr << "usage: " << std::string(argv[0])
                << " plot-time final-time dt objects" << std::endl
                << " Details:" << std::endl
                << " ----------------------------------" << std::endl
                << "  plot-time:       interval after how many time units to plot."
                << " Use 0 to switch off plotting" << std::endl
                << "  final-time:      simulated time (greater 0)" << std::endl
                << "  dt:              time step size (greater 0)" << std::endl
                << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
                << std::endl;
      throw -1;
    }
    else if ((argc-4)%7!=0) {
      std::cerr << "error in arguments: each body is given by seven entries"
                << " (position, velocity, mass)" << std::endl;
      std::cerr << "got " << argc << " arguments"
                << " (three of them are reserved)" << std::endl;
      throw -2;
    }
  }

  // Helper function to detect and merge colliding bodies
  void mergeBodies() {
    bool mergeOccurred = false;

    for (int i = 0; i < NumberOfBodies; i++) {
      if (merged[i]) continue;  // Skip already merged bodies

      for (int j = i + 1; j < NumberOfBodies; j++) {
        if (merged[j]) continue;  // Skip already merged bodies

        // Calculate the Euclidean distance between bodies i and j
        const double dx = bodies[j].x - bodies[i].x;
        const double dy = bodies[j].y - bodies[i].y;
        const double dz = bodies[j].z - bodies[i].z;
        const double distance = sqrt(dx*dx + dy*dy + dz*dz);

        // Check if they should merge: |x1-x2| ≤ C(m1+m2)
        if (distance <= collisionConstant * (bodies[i].mass + bodies[j].mass)) {
          mergeOccurred = true;

          // Calculate the new position and velocity as the mass-weighted mean
          double totalMass = bodies[i].mass + bodies[j].mass;
          
          // Position: x = (m1*x1 + m2*x2)/(m1+m2)
          bodies[i].x = (bodies[i].mass * bodies[i].x + bodies[j].mass * bodies[j].x) / totalMass;
          bodies[i].y = (bodies[i].mass * bodies[i].y + bodies[j].mass * bodies[j].y) / totalMass;
          bodies[i].z = (bodies[i].mass * bodies[i].z + bodies[j].mass * bodies[j].z) / totalMass;
          
          // Velocity: v = (m1*v1 + m2*v2)/(m1+m2)
          bodies[i].vx = (bodies[i].mass * bodies[i].vx + bodies[j].mass * bodies[j].vx) / totalMass;
          bodies[i].vy = (bodies[i].mass * bodies[i].vy + bodies[j].mass * bodies[j].vy) / totalMass;
          bodies[i].vz = (bodies[i].mass * bodies[i].vz + bodies[j].mass * bodies[j].vz) / totalMass;

          // Update the mass of body i
          bodies[i].mass = totalMass;

          // Mark body j as merged
          merged[j] = true;

          std::cout << "Bodies " << i << " and " << j << " merged with new mass "
                    << bodies[i].mass << " at position ("
                    << bodies[i].x << ", " << bodies[i].y << ", " << bodies[i].z << ")" << std::endl;
        }
      }
    }

    // If merges occurred, we need to compact the bodies array
    if (mergeOccurred) {
      int newIndex = 0;
      for (int i = 0; i < NumberOfBodies; i++) {
        if (!merged[i]) {
          // If this body wasn't merged, move it to the next available slot
          if (newIndex != i) {
            bodies[newIndex] = bodies[i];
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
 */
class NBodyBarnsHutSimulation : public NBodySimulationCollision {
private:
    OctreeNode* octree;       // The octree for Barnes-Hut algorithm
    double domainSize;        // Size of the simulation domain
    
    // Adapter arrays for compatibility with OctreeNode
    double** positionAdapter;
    double* massAdapter;
    
public:
    // Constructor
    NBodyBarnsHutSimulation() : 
        NBodySimulationCollision(), 
        octree(nullptr), 
        domainSize(0.0),
        positionAdapter(nullptr),
        massAdapter(nullptr) {}
    
    // Destructor
    ~NBodyBarnsHutSimulation() {
        if (octree != nullptr) {
            delete octree;
            octree = nullptr;
        }
        
        // Clean up adapter arrays
        cleanupAdapters();
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
    void updateBody() override {
        try {
            // Initial setup
            timeStepCounter++;
            maxV = 0.0;
            minDx = std::numeric_limits<double>::max();

            // Create adapter arrays for OctreeNode
            createAdapters();
            
            // Build the octree for this time step
            buildOctree();

            // Allocate and initialize force arrays
            double* force0 = new double[NumberOfBodies]();
            double* force1 = new double[NumberOfBodies]();
            double* force2 = new double[NumberOfBodies]();

            // Compute forces for each body using Barnes-Hut algorithm
            for (int i = 0; i < NumberOfBodies; i++) {
                // Skip bodies with invalid positions
                if (std::isnan(bodies[i].x) || std::isinf(bodies[i].x) ||
                    std::isnan(bodies[i].y) || std::isinf(bodies[i].y) ||
                    std::isnan(bodies[i].z) || std::isinf(bodies[i].z)) {
                    continue;
                }
                
                // Temporary array for forces from octree
                double forces[3] = {0.0, 0.0, 0.0};
                
                // Compute force using the octree
                if (octree != nullptr) {
                    octree->computeForce(i, positionAdapter, forces, massAdapter);
                }
                
                // Store the forces
                force0[i] = forces[0];
                force1[i] = forces[1];
                force2[i] = forces[2];
            }

            // Update positions using current velocities
            for (int i = 0; i < NumberOfBodies; i++) {
                bodies[i].x += timeStepSize * bodies[i].vx;
                bodies[i].y += timeStepSize * bodies[i].vy;
                bodies[i].z += timeStepSize * bodies[i].vz;
                
                // Update adapters after position changes
                updateAdapters();
            }

            // Check for collisions and merge bodies if needed
            mergeBodies();
            
            // Update adapters after possible merges
            updateAdapters();

            // Update velocities using computed forces
            for (int i = 0; i < NumberOfBodies; i++) {
                bodies[i].vx += timeStepSize * force0[i] / bodies[i].mass;
                bodies[i].vy += timeStepSize * force1[i] / bodies[i].mass;
                bodies[i].vz += timeStepSize * force2[i] / bodies[i].mass;

                // Update maximum velocity
                maxV = std::max(maxV, std::sqrt(bodies[i].vx*bodies[i].vx + 
                                              bodies[i].vy*bodies[i].vy + 
                                              bodies[i].vz*bodies[i].vz));
            }

            // Recalculate minimum distance between bodies
            if (NumberOfBodies > 1) {
                for (int i = 0; i < NumberOfBodies - 1; i++) {
                    for (int j = i + 1; j < NumberOfBodies; j++) {
                        double dx = bodies[j].x - bodies[i].x;
                        double dy = bodies[j].y - bodies[i].y;
                        double dz = bodies[j].z - bodies[i].z;
                        double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
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
    // Create adapter arrays for compatibility with OctreeNode
    void createAdapters() {
        // Clean up old adapters if they exist
        cleanupAdapters();
        
        // Create new adapters
        positionAdapter = new double*[NumberOfBodies];
        massAdapter = new double[NumberOfBodies];
        
        for (int i = 0; i < NumberOfBodies; i++) {
            positionAdapter[i] = new double[3];
            positionAdapter[i][0] = bodies[i].x;
            positionAdapter[i][1] = bodies[i].y;
            positionAdapter[i][2] = bodies[i].z;
            massAdapter[i] = bodies[i].mass;
        }
    }
    
    // Update adapters to match the current state of bodies
    void updateAdapters() {
        if (positionAdapter == nullptr || massAdapter == nullptr) {
            createAdapters();
            return;
        }
        
        for (int i = 0; i < NumberOfBodies; i++) {
            positionAdapter[i][0] = bodies[i].x;
            positionAdapter[i][1] = bodies[i].y;
            positionAdapter[i][2] = bodies[i].z;
            massAdapter[i] = bodies[i].mass;
        }
    }
    
    // Clean up adapter arrays
    void cleanupAdapters() {
        if (positionAdapter != nullptr) {
            for (int i = 0; i < NumberOfBodies; i++) {
                delete[] positionAdapter[i];
            }
            delete[] positionAdapter;
            positionAdapter = nullptr;
        }
        
        if (massAdapter != nullptr) {
            delete[] massAdapter;
            massAdapter = nullptr;
        }
    }
    
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
                // Skip bodies with invalid positions
                if (std::isnan(bodies[i].x) || std::isinf(bodies[i].x) ||
                    std::isnan(bodies[i].y) || std::isinf(bodies[i].y) ||
                    std::isnan(bodies[i].z) || std::isinf(bodies[i].z)) {
                    continue;
                }
                
                // Insert into the octree using adapter arrays
                octree->insert(i, positionAdapter, massAdapter);
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
                if (!std::isnan(bodies[i].x) && !std::isinf(bodies[i].x))
                    maxCoord = std::max(maxCoord, std::abs(bodies[i].x));
                
                if (!std::isnan(bodies[i].y) && !std::isinf(bodies[i].y))
                    maxCoord = std::max(maxCoord, std::abs(bodies[i].y));
                
                if (!std::isnan(bodies[i].z) && !std::isinf(bodies[i].z))
                    maxCoord = std::max(maxCoord, std::abs(bodies[i].z));
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