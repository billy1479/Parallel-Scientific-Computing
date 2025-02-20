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
  
  std::cout << "Checking time step quality..." << std::endl;
  checkTimeStepQuality();
  std::cout << "Finding optimal time step..." << std::endl;
  
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
    
    // // Calculate and store time step metrics
    // double energyAfter = calculateTotalEnergy();
    // double energyError = std::abs((energyAfter - energyBefore) / energyBefore);
    // double relativeError = std::abs((energyAfter - energyBefore) / energyBefore);
    // timeStepHistory.push_back({timeStepSize, energyError, t});

    // if(relativeError > 1e-3) {  // More than 0.1% change
    //   std::cout << "Warning: Large energy change detected at step " 
    //             << timeStepCounter << ": " << relativeError * 100 
    //             << "%" << std::endl;
    // }
    
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

void NBodySimulation::calculateOptimalTimeStep() {
  // Find shortest orbital period among all pairs
  double minPeriod = std::numeric_limits<double>::max();
  for(int i = 0; i < NumberOfBodies; i++) {
      for(int j = i+1; j < NumberOfBodies; j++) {
          double period = calculateOrbitalPeriod(i, j);
          minPeriod = std::min(minPeriod, period);
      }
  }
  
  // Calculate system's dynamical time
  double tDynamical = calculateDynamicalTime();
  
  // The optimal time step should be small enough to:
  // 1. Resolve the shortest orbital period (need ~100 steps per orbit)
  // 2. Resolve the system's dynamical time (need ~100 steps per crossing)
  double dtFromOrbits = minPeriod / 100.0;
  double dtFromDynamics = tDynamical / 100.0;
  
  // Choose the more restrictive condition
  double recommendedTimeStep = std::min(dtFromOrbits, dtFromDynamics);
  
  std::cout << "\nTime Step Analysis Based on Physical Timescales:" << std::endl;
  std::cout << "=============================================" << std::endl;
  std::cout << "Shortest orbital period: " << minPeriod << std::endl;
  std::cout << "System dynamical time: " << tDynamical << std::endl;
  std::cout << "Recommended time step: " << recommendedTimeStep << std::endl;
  std::cout << "\nThis recommendation is based on:" << std::endl;
  std::cout << "- Resolving the fastest orbital motion in the system" << std::endl;
  std::cout << "- Ensuring accurate integration of system-wide dynamics" << std::endl;
  
  // Provide context for the recommendation
  if(timeStepSize > recommendedTimeStep) {
      double ratio = timeStepSize / recommendedTimeStep;
      std::cout << "\nCurrent time step is " << ratio 
                << " times larger than recommended." << std::endl;
      std::cout << "This explains the observed energy conservation issues." << std::endl;
  } else {
      std::cout << "\nCurrent time step is within acceptable range." << std::endl;
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
            
  // if (!timeStepHistory.empty()) {
  //     auto bestMetrics = std::min_element(
  //         timeStepHistory.begin(),
  //         timeStepHistory.end(),
  //         [](const TimeStepMetrics& a, const TimeStepMetrics& b) -> bool {
  //             return a.energyError < b.energyError;
  //         }
  //     );
      
  //     // Calculate average energy error
  //     double avgError = 0.0;
  //     for (const auto& metrics : timeStepHistory) {
  //         avgError += metrics.energyError;
  //     }
  //     avgError /= timeStepHistory.size();
      
  //     std::cout << "\nTime Step Analysis:" << std::endl;
  //     std::cout << "-------------------" << std::endl;
  //     std::cout << "Initial time step: " << timeStepHistory.front().timeStep << std::endl;
  //     std::cout << "Best time step found: " << bestMetrics->timeStep << std::endl;
  //     std::cout << "Best energy error: " << bestMetrics->energyError << std::endl;
  //     std::cout << "Average energy error: " << avgError << std::endl;
  //     std::cout << "Time of best performance: " << bestMetrics->time << std::endl;
      
  //     // Provide recommendations based on energy error
  //     std::cout << "\nRecommendation:" << std::endl;
  //     if (bestMetrics->energyError < 1e-6) {
  //         std::cout << "Time step performance is excellent - current time step is well-suited for this simulation" << std::endl;
  //     } else if (bestMetrics->energyError < 1e-4) {
  //         std::cout << "Time step performance is good - current time step provides good accuracy" << std::endl;
  //     } else {
  //         std::cout << "Consider using a smaller time step (try " 
  //                   << bestMetrics->timeStep * 0.5 
  //                   << ") for better accuracy" << std::endl;
  //     }

  //     // Add energy conservation analysis
  //   double finalEnergy = calculateTotalEnergy();
  //   double totalEnergyError = std::abs((finalEnergy - initialEnergy) / initialEnergy);
    
  //   std::cout << "\nEnergy Conservation Analysis:" << std::endl;
  //   std::cout << "============================" << std::endl;
  //   std::cout << "Initial Energy: " << initialEnergy << std::endl;
  //   std::cout << "Final Energy: " << finalEnergy << std::endl;
  //   std::cout << "Total Energy Error: " << totalEnergyError * 100 << "%" << std::endl;
    
  //   // Provide interpretation of results
  //   if(totalEnergyError < 1e-3) {
  //       std::cout << "\nTime Step Assessment: EXCELLENT" << std::endl;
  //       std::cout << "The current time step provides very good energy conservation." << std::endl;
  //   } else if(totalEnergyError < 1e-2) {
  //       std::cout << "\nTime Step Assessment: ACCEPTABLE" << std::endl;
  //       std::cout << "The time step gives reasonable results but could be reduced for better accuracy." << std::endl;
  //   } else {
  //       std::cout << "\nTime Step Assessment: NEEDS ADJUSTMENT" << std::endl;
  //       std::cout << "Consider reducing the time step to improve accuracy." << std::endl;
  //       std::cout << "Suggested time step: " << timeStepSize * 0.5 << std::endl;
  //   }
  // }
}