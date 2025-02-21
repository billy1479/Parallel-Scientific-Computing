#include "NBodySimulation.h"
#include <algorithm>

NBodySimulation::NBodySimulation() :
    t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
    x(nullptr), v(nullptr), mass(nullptr),
    timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
    snapshotCounter(0), timeStepCounter(0), initialEnergy(0) {
}

NBodySimulation::~NBodySimulation () {
  if (x != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] x[i];
    delete [] x;
  }
  if (v != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] v[i];
    delete [] v;
  }
  if (mass != nullptr) {
    delete [] mass;
  }
}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp(int argc, char** argv) {
  checkInput(argc, argv);
  
  NumberOfBodies = (argc-4) / 7;
  
  // Allocate arrays
  x = new double*[NumberOfBodies];
  v = new double*[NumberOfBodies];
  mass = new double[NumberOfBodies];
  
  int readArgument = 1;
  
  tPlotDelta = std::stof(argv[readArgument++]);
  tFinal = std::stof(argv[readArgument++]);
  timeStepSize = std::stof(argv[readArgument++]);

  std::cout << "tPlotDelta: " << tPlotDelta << std::endl;
  std::cout << "tFinal: " << tFinal << std::endl;
  std::cout << "timeStepSize: " << timeStepSize << std::endl;
  
  // Initialize body data
  for (int i=0; i<NumberOfBodies; i++) {
      x[i] = new double[3];
      v[i] = new double[3];
      
      // Read position
      x[i][0] = std::stof(argv[readArgument++]);
      x[i][1] = std::stof(argv[readArgument++]);
      x[i][2] = std::stof(argv[readArgument++]);
      
      // Read velocity
      v[i][0] = std::stof(argv[readArgument++]);
      v[i][1] = std::stof(argv[readArgument++]);
      v[i][2] = std::stof(argv[readArgument++]);
      
      // Read mass
      mass[i] = std::stof(argv[readArgument++]);
      
      if (mass[i] <= 0.0) {
          std::cerr << "invalid mass for body " << i << std::endl;
          exit(-2);
      }
  }
  std::cout << "Calculate time step and assign value..." << std::endl;
  calculateComprehensiveTimeStep();
  std::cout << "Checking time step quality..." << std::endl;
  checkTimeStepQuality();
  // std::cout << "Finding optimal time step..." << std::endl;
  std::cout << "Measuring convergence order..." << std::endl;
  measureConvergenceOrder();
  
  // Configure plotting
  if (tPlotDelta <= 0.0) {
      std::cout << "plotting switched off" << std::endl;
      tPlot = tFinal + 1.0;
  } else {
      std::cout << "plot initial setup plus every " << tPlotDelta 
                << " time units" << std::endl;
      tPlot = 0.0;
  }
}

void NBodySimulation::adaptiveTimeStep() {
  // CFL-like condition
  double safetyFactor = 0.1;  // Conservative starting point
  double newTimeStep = safetyFactor * minDx / maxV;
  
  // Don't change time step too abruptly
  double maxChange = 1.5;  // Allow 50% increase at most
  timeStepSize = std::min(
      newTimeStep,
      timeStepSize * maxChange
  );
}

// Add this function to your simulation class
void NBodySimulation::findOptimalTimeStep(double initialGuess) {
  double dt = initialGuess;
  double prevEnergy = calculateTotalEnergy();
  
  for(int i = 0; i < 5; i++) {  // Test different time steps
      timeStepSize = dt;
      double newEnergy = calculateTotalEnergy();
      double relativeError = std::abs((newEnergy - prevEnergy)/prevEnergy);
      
      std::cout << "dt: " << dt << " Error: " << relativeError << std::endl;
      
      dt /= 2.0;  // Halve the time step each iteration
  }
}

double NBodySimulation::calculateTotalEnergy() {
  double energy = 0.0;
  // Kinetic energy
  for(int i = 0; i < NumberOfBodies; i++) {
      double velocitySquared = v[i][0]*v[i][0] + 
                             v[i][1]*v[i][1] + 
                             v[i][2]*v[i][2];
      energy += 0.5 * mass[i] * velocitySquared;
  }
  // Potential energy
  for(int i = 0; i < NumberOfBodies; i++) {
      for(int j = i+1; j < NumberOfBodies; j++) {
          double distance = sqrt(
              (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
              (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
              (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
          );
          energy -= (mass[i] * mass[j]) / distance;
      }
  }
  return energy;
}

void NBodySimulation::printEnergySummary() {
  double currentEnergy = calculateTotalEnergy();
  double relativeError = std::abs((currentEnergy - initialEnergy) / initialEnergy);
  
  std::cout << "Energy Analysis:" << std::endl;
  std::cout << "Initial Energy: " << initialEnergy << std::endl;
  std::cout << "Current Energy: " << currentEnergy << std::endl;
  std::cout << "Relative Error: " << relativeError << std::endl;
  
  // For a symplectic integrator like yours, the energy error should:
  // 1. Remain bounded (not grow exponentially)
  // 2. Typically stay below 10^-3 for a good time step
  if(relativeError < 1e-3) {
      std::cout << "Time step appears appropriate (energy error < 0.1%)" << std::endl;
  } else if(relativeError < 1e-2) {
      std::cout << "Time step may need reduction (energy error < 1%)" << std::endl;
  } else {
      std::cout << "Time step likely too large (energy error > 1%)" << std::endl;
  }
}

double NBodySimulation::force_calculation (int i, int j, int direction){
  // Euclidean distance
  const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );
                               
  const double distance3 = distance * distance * distance;
  minDx = std::min( minDx,distance );

  double temp = (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;

  return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
}

void NBodySimulation::updateBody() {
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();
    
    // Store energy before update
    // double energyBefore = calculateTotalEnergy();
    
    // Allocate and initialize force arrays
    double* force0 = new double[NumberOfBodies]();
    double* force1 = new double[NumberOfBodies]();
    double* force2 = new double[NumberOfBodies]();
    
    // Special case for single body
    if (NumberOfBodies == 1) minDx = 0;
    
    // Calculate forces between all pairs of bodies
    for (int i=0; i<NumberOfBodies; i++) {
        for (int j=i+1; j<NumberOfBodies; j++) {
            if (i != j) {
                force0[i] += force_calculation(i,j,0);
                force1[i] += force_calculation(i,j,1);
                force2[i] += force_calculation(i,j,2);
                
                force0[j] -= force_calculation(i,j,0);
                force1[j] -= force_calculation(i,j,1);
                force2[j] -= force_calculation(i,j,2);
            }
        }
    }
    
    // Update positions (first part of symplectic Euler)
    for (int i=0; i<NumberOfBodies; i++) {
        x[i][0] = x[i][0] + timeStepSize * v[i][0];
        x[i][1] = x[i][1] + timeStepSize * v[i][1];
        x[i][2] = x[i][2] + timeStepSize * v[i][2];
    }
    
    // Update velocities (second part of symplectic Euler)
    for (int i=0; i<NumberOfBodies; i++) {
        v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
        v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
        v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];
        
        maxV = std::max(maxV, std::sqrt(
            v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]
        ));
    }
  
    // Update simulation time
    t += timeStepSize;
    
    // Clean up
    delete[] force0;
    delete[] force1;
    delete[] force2;
}

bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}

void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

double NBodySimulation::calculateOrbitalPeriod(int i, int j) {
  // Calculate distance between bodies i and j
  double distance = sqrt(
      (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
      (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
      (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
  );
  
  // Calculate relative velocity
  double relativeV = sqrt(
      pow(v[j][0]-v[i][0], 2) +
      pow(v[j][1]-v[i][1], 2) +
      pow(v[j][2]-v[i][2], 2)
  );
  
  // For a bound pair, the orbital period depends on separation and masses
  // T = 2π * sqrt(r³/(G*(m1+m2)))
  // Note: In gravitational units, G = 1
  return 2.0 * M_PI * sqrt(
      pow(distance, 3) / (mass[i] + mass[j])
  );
}

double NBodySimulation::calculateDynamicalTime() {
  // The dynamical time is characteristic of the whole system
  // It's approximately the time for a typical particle to cross the system
  
  // First find system size
  double maxR = 0.0;
  double totalMass = 0.0;
  
  for(int i = 0; i < NumberOfBodies; i++) {
      double r = sqrt(
          x[i][0]*x[i][0] + 
          x[i][1]*x[i][1] + 
          x[i][2]*x[i][2]
      );
      maxR = std::max(maxR, r);
      totalMass += mass[i];
  }
  
  // Dynamical time = sqrt(R³/(GM))
  return sqrt(pow(maxR, 3) / totalMass);
}

void NBodySimulation::calculateComprehensiveTimeStep() {
    // This method combines multiple approaches to determine the optimal time step
    // by considering initial conditions, system properties, and numerical stability
    
    // ===== PART 1: PHYSICAL TIMESCALES =====
    
    // Calculate characteristic timescales of the system
    double minOrbitalPeriod = std::numeric_limits<double>::max();
    double totalMass = 0.0;
    double maxDistance = 0.0;
    double minDistance = std::numeric_limits<double>::max();
    double maxVelocity = 0.0;
    double maxAcceleration = 0.0;
    
    // Calculate center of mass
    double comX = 0.0, comY = 0.0, comZ = 0.0;
    for (int i = 0; i < NumberOfBodies; i++) {
        totalMass += mass[i];
        comX += mass[i] * x[i][0];
        comY += mass[i] * x[i][1];
        comZ += mass[i] * x[i][2];
    }
    comX /= totalMass;
    comY /= totalMass;
    comZ /= totalMass;
    
    // Analyze all body pairs for orbital periods and distances
    for (int i = 0; i < NumberOfBodies; i++) {
        // Distance from center of mass (for system size estimate)
        double distFromCOM = sqrt(
            pow(x[i][0] - comX, 2) +
            pow(x[i][1] - comY, 2) +
            pow(x[i][2] - comZ, 2)
        );
        maxDistance = std::max(maxDistance, distFromCOM);
        
        // Current velocity magnitude
        double vel = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
        maxVelocity = std::max(maxVelocity, vel);
        
        // Calculate max acceleration from all pairwise forces
        for (int j = 0; j < NumberOfBodies; j++) {
            if (i == j) continue;
            
            // Calculate distance between bodies
            double dx = x[j][0] - x[i][0];
            double dy = x[j][1] - x[i][1];
            double dz = x[j][2] - x[i][2];
            double distSquared = dx*dx + dy*dy + dz*dz;
            double dist = sqrt(distSquared);
            
            minDistance = std::min(minDistance, dist);
            
            // For close pairs, calculate orbital period
            if (i < j) {  // avoid double counting
                // Two-body orbital period (Kepler's third law)
                double period = 2.0 * M_PI * sqrt(
                    pow(dist, 3) / (mass[i] + mass[j])
                );
                minOrbitalPeriod = std::min(minOrbitalPeriod, period);
            }
            
            // Calculate acceleration magnitude due to body j
            double acc = mass[j] / (distSquared * dist);
            maxAcceleration = std::max(maxAcceleration, acc);
        }
    }
    
    // ===== PART 2: NUMERICAL STABILITY CRITERIA =====
    
    // Safety factor (more conservative for systems with many bodies)
    double safetyFactor = 0.05 * pow(0.99, std::min(NumberOfBodies / 100, 10)); 
    
    // 1. CFL-like condition: dt < dx/v
    double dtCFL = std::numeric_limits<double>::max();
    if (maxVelocity > 1e-10) {  // Prevent division by zero
        dtCFL = safetyFactor * minDistance / maxVelocity;
    }
    
    // 2. Acceleration-based criterion: dt² < dx/a
    double dtAccel = std::numeric_limits<double>::max();
    if (maxAcceleration > 1e-10) {  // Prevent division by zero
        dtAccel = sqrt(safetyFactor * minDistance / maxAcceleration);
    }
    
    // 3. Orbital resolution criterion: need ~100 steps per orbit
    double dtOrbital = minOrbitalPeriod / 100.0;
    
    // 4. System dynamical time (crossing time)
    double dynTime = sqrt(pow(maxDistance, 3) / totalMass);
    double dtDynamical = dynTime / 100.0;
    
    // 5. Energy conservation - use fewer test steps to avoid memory issues
    double initialEnergy = calculateTotalEnergy();
    double dtEnergy = std::min(dtAccel, dtOrbital);  // Start with a physics-based estimate
    
    // Only test 2-3 carefully chosen timesteps based on the physics criteria
    std::vector<double> testSteps;
    
    // Add the smallest physics-based timestep
    double smallestPhysicsStep = std::min({
        dtAccel, 
        dtOrbital,
        dtDynamical,
        tFinal / 1000.0  // Cap very small steps
    });
    
    // Make sure our test steps are valid
    if (smallestPhysicsStep <= 0 || !std::isfinite(smallestPhysicsStep)) {
        smallestPhysicsStep = tFinal / 1000.0;
    }
    
    // Add only a few carefully chosen test steps
    testSteps.push_back(smallestPhysicsStep);
    testSteps.push_back(smallestPhysicsStep * 2.0);
    
    // For safety with large systems (over 100 bodies), use a simpler approach
    if (NumberOfBodies > 100) {
        // For large systems, we'll skip the energy testing and use physics-based criteria only
        dtEnergy = smallestPhysicsStep;
        std::cout << "Large system detected: Using physics-based time step without energy testing" << std::endl;
    } else {
        // For smaller systems, we can do the full test
        // Save initial state using stack-based storage for small arrays
        std::vector<std::vector<std::vector<double>>> xBackup(NumberOfBodies, 
            std::vector<std::vector<double>>(1, std::vector<double>(3)));
        std::vector<std::vector<std::vector<double>>> vBackup(NumberOfBodies, 
            std::vector<std::vector<double>>(1, std::vector<double>(3)));
        
        // Save current state
        for (int i = 0; i < NumberOfBodies; i++) {
            for (int d = 0; d < 3; d++) {
                xBackup[i][0][d] = x[i][d];
                vBackup[i][0][d] = v[i][d];
            }
        }
        
        double bestError = std::numeric_limits<double>::max();
        
        // Test each timestep
        for (double testStep : testSteps) {
            // Restore initial state
            for (int i = 0; i < NumberOfBodies; i++) {
                for (int d = 0; d < 3; d++) {
                    x[i][d] = xBackup[i][0][d];
                    v[i][d] = vBackup[i][0][d];
                }
            }
            
            // Test this timestep for just a few steps
            double originalT = t;
            timeStepSize = testStep;
            
            // For smaller systems, we can afford more test steps
            int numTestSteps = std::min(5, 1000 / NumberOfBodies);
            for (int step = 0; step < numTestSteps; step++) {
                updateBody();
            }
            
            // Evaluate energy conservation
            double finalEnergy = calculateTotalEnergy();
            double error = std::abs((finalEnergy - initialEnergy) / initialEnergy);
            
            std::cout << "Test dt: " << testStep << " -> Energy error: " << error << std::endl;
            
            if (error < bestError) {
                bestError = error;
                dtEnergy = testStep;
            }
            
            // Reset time
            t = originalT;
        }
        
        // Restore initial state one last time
        for (int i = 0; i < NumberOfBodies; i++) {
            for (int d = 0; d < 3; d++) {
                x[i][d] = xBackup[i][0][d];
                v[i][d] = vBackup[i][0][d];
            }
        }
    }
    
    // ===== PART 3: DETERMINE FINAL TIME STEP =====
    
    // Validate all time steps (replace infinities, NaNs, and zeros with safe values)
    if (!std::isfinite(dtCFL) || dtCFL <= 0) dtCFL = tFinal / 100.0;
    if (!std::isfinite(dtAccel) || dtAccel <= 0) dtAccel = tFinal / 100.0;
    if (!std::isfinite(dtOrbital) || dtOrbital <= 0) dtOrbital = tFinal / 100.0;
    if (!std::isfinite(dtDynamical) || dtDynamical <= 0) dtDynamical = tFinal / 100.0;
    if (!std::isfinite(dtEnergy) || dtEnergy <= 0) dtEnergy = tFinal / 100.0;
    
    // Gather all criteria
    std::vector<std::pair<double, std::string>> criteria = {
        {dtCFL, "CFL condition"},
        {dtAccel, "Acceleration limit"},
        {dtOrbital, "Shortest orbital period / 100"},
        {dtDynamical, "System dynamical time / 100"},
        {dtEnergy, "Empirical energy conservation"}
    };
    
    // Find the limiting criterion (smallest time step)
    double finalTimeStep = tFinal / 10.0;  // Start with a reasonable fraction of total time
    std::string limitingFactor = "Default (tFinal/10)";
    
    for (const auto& criterion : criteria) {
        if (criterion.first > 0 && criterion.first < finalTimeStep) {
            finalTimeStep = criterion.first;
            limitingFactor = criterion.second;
        }
    }
    
    // Additional safety for large systems
    if (NumberOfBodies > 500) {
        double safetyScale = 0.8;  // 20% smaller for large systems
        finalTimeStep *= safetyScale;
        std::cout << "Applied additional " << (1-safetyScale)*100 << "% safety factor for large system" << std::endl;
    }
    
    // Ensure time step isn't too small compared to simulation length
    if (finalTimeStep < tFinal / 1e6) {
        std::cout << "Warning: Calculated time step is extremely small." << std::endl;
        std::cout << "Consider modifying initial conditions or using a specialized integrator." << std::endl;
        finalTimeStep = tFinal / 1e6;  // Cap at a reasonable minimum
    }
    
    // Ensure time step isn't too large compared to total simulation time
    if (finalTimeStep > tFinal / 10.0) {
        finalTimeStep = tFinal / 10.0;
        std::cout << "Warning: Calculated time step is very large. Limiting to tFinal/10." << std::endl;
        limitingFactor = "Maximum limit (tFinal/10)";
    }
    
    // ===== PART 4: OUTPUT RESULTS =====
    
    std::cout << "\n==== COMPREHENSIVE TIME STEP ANALYSIS ====" << std::endl;
    std::cout << "System properties:" << std::endl;
    std::cout << "  * Bodies: " << NumberOfBodies << std::endl;
    std::cout << "  * Total mass: " << totalMass << std::endl;
    std::cout << "  * Max system extent: " << maxDistance << std::endl;
    std::cout << "  * Min separation: " << minDistance << std::endl;
    std::cout << "  * Max velocity: " << maxVelocity << std::endl;
    std::cout << "  * Shortest orbital period: " << minOrbitalPeriod << std::endl;
    std::cout << "  * System dynamical time: " << dynTime << std::endl;
    
    std::cout << "\nTime step criteria:" << std::endl;
    for (const auto& criterion : criteria) {
        std::cout << "  * " << criterion.second << ": " << criterion.first;
        if (criterion.second == limitingFactor) {
            std::cout << " (LIMITING FACTOR)";
        }
        std::cout << std::endl;
    }
    
    std::cout << "\nFinal time step recommendation: " << finalTimeStep << std::endl;
    std::cout << "Based on: " << limitingFactor << std::endl;
    
    // Return information to caller and update timeStepSize
    timeStepSize = finalTimeStep;
    
    // Calculate and store initial energy for future reference
    initialEnergy = calculateTotalEnergy();
    
    // Estimate computational requirements
    long estimatedSteps = static_cast<long>(tFinal / timeStepSize);
    std::cout << "\nEstimated number of time steps: " << estimatedSteps << std::endl;
    
    // Provide warnings if the computation might be too intensive
    if (estimatedSteps > 1000000) {
        std::cout << "WARNING: This simulation will require a very large number of steps." << std::endl;
        std::cout << "Consider using a larger time step or shorter simulation time." << std::endl;
    }
    
    if (tPlotDelta > 0) {
        int plotSteps = int(tFinal / tPlotDelta);
        std::cout << "Will generate approximately " << plotSteps << " snapshots" << std::endl;
        
        // Warn about excessive file output
        if (plotSteps > 1000) {
            std::cout << "WARNING: Large number of output files will be generated." << std::endl;
            std::cout << "Consider increasing tPlotDelta." << std::endl;
        }
    }
    
    // Estimate memory usage (rough approximation)
    double memoryPerBody = 6 * sizeof(double) + sizeof(double*) * 2;  // positions, velocities, mass
    double totalMemoryMB = (NumberOfBodies * memoryPerBody) / (1024 * 1024);
    std::cout << "Estimated memory usage: " << totalMemoryMB << " MB" << std::endl;
}

void NBodySimulation::measureConvergenceOrder() {
  // Store the original time step
  double originalTimeStep = timeStepSize;
  
  // Set up a reference solution with a very small time step
  double referenceTimeStep = originalTimeStep / 64.0;
  
  // Save initial state
  std::vector<std::vector<double>> initialX(NumberOfBodies, std::vector<double>(3));
  std::vector<std::vector<double>> initialV(NumberOfBodies, std::vector<double>(3));
  
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          initialX[i][d] = x[i][d];
          initialV[i][d] = v[i][d];
      }
  }
  
  // Time to simulate for convergence testing
  double testTime = 1.0;  // Simulate for 1 time unit
  
  // Reference solution
  timeStepSize = referenceTimeStep;
  int refSteps = static_cast<int>(testTime / timeStepSize);
  
  for (int step = 0; step < refSteps; step++) {
      updateBody();
  }
  
  // Store reference solution
  std::vector<std::vector<double>> referenceX(NumberOfBodies, std::vector<double>(3));
  
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          referenceX[i][d] = x[i][d];
      }
  }
  
  // Test with different time steps
  std::vector<double> timeSteps;
  std::vector<double> errors;
  
  // Use time steps that are powers of 2 times the reference
  timeSteps.push_back(referenceTimeStep * 2);   // h/32
  timeSteps.push_back(referenceTimeStep * 4);   // h/16
  timeSteps.push_back(referenceTimeStep * 8);   // h/8
  timeSteps.push_back(referenceTimeStep * 16);  // h/4
  timeSteps.push_back(referenceTimeStep * 32);  // h/2
  
  for (double testTimeStep : timeSteps) {
      // Reset to initial state
      for (int i = 0; i < NumberOfBodies; i++) {
          for (int d = 0; d < 3; d++) {
              x[i][d] = initialX[i][d];
              v[i][d] = initialV[i][d];
          }
      }
      
      // Run with current time step
      timeStepSize = testTimeStep;
      int numSteps = static_cast<int>(testTime / timeStepSize);
      
      for (int step = 0; step < numSteps; step++) {
          updateBody();
      }
      
      // Measure error (RMS of position difference)
      double sumSquaredError = 0.0;
      for (int i = 0; i < NumberOfBodies; i++) {
          for (int d = 0; d < 3; d++) {
              double err = x[i][d] - referenceX[i][d];
              sumSquaredError += err * err;
          }
      }
      double rmsError = sqrt(sumSquaredError / (NumberOfBodies * 3));
      errors.push_back(rmsError);
      
      std::cout << "Time step: " << testTimeStep 
                << ", Error: " << rmsError << std::endl;
  }
  
  // Calculate convergence order using log-log linear regression
  double sumLogTimeStep = 0.0;
  double sumLogError = 0.0;
  double sumLogTimeStepSquared = 0.0;
  double sumLogTimeStepLogError = 0.0;
  int n = timeSteps.size();
  
  for (int i = 0; i < n; i++) {
      double logTimeStep = log(timeSteps[i]);
      double logError = log(errors[i]);
      
      sumLogTimeStep += logTimeStep;
      sumLogError += logError;
      sumLogTimeStepSquared += logTimeStep * logTimeStep;
      sumLogTimeStepLogError += logTimeStep * logError;
  }
  
  // Slope of the best-fit line gives the convergence order
  double convergenceOrder = (n * sumLogTimeStepLogError - sumLogTimeStep * sumLogError) /
                            (n * sumLogTimeStepSquared - sumLogTimeStep * sumLogTimeStep);
  
  std::cout << "Measured convergence order: " << convergenceOrder << std::endl;
  
  // Restore original time step
  timeStepSize = originalTimeStep;
  
  // Reset to initial state
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          x[i][d] = initialX[i][d];
          v[i][d] = initialV[i][d];
      }
  }
}

void NBodySimulation::checkTimeStepQuality() {
  // Track energy changes over time
  double initialEnergy = calculateTotalEnergy();
  double prevEnergy = initialEnergy;
  bool hasLargeFluctuations = false;
  
  // Run simulation for a short test period
  for(int i = 0; i < 100; i++) {  // Test over 100 steps
      updateBody();
      double currentEnergy = calculateTotalEnergy();
      double relativeChange = std::abs((currentEnergy - prevEnergy)/prevEnergy);
      
      if(relativeChange > 0.01) {  // 1% threshold
          std::cout << "Warning: Large energy change detected at step " << i << ": " << relativeChange * 100 << "%" << std::endl;
          hasLargeFluctuations = true;
          break;
      }
      prevEnergy = currentEnergy;
  }
}

void NBodySimulation::printSummary() {
  std::cout << "Simulation Summary:" << std::endl;
  std::cout << "==================" << std::endl;
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
            
}