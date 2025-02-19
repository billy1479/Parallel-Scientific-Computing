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
  
  // Calculate and store initial energy
  initialEnergy = calculateTotalEnergy();
  
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
    double energyBefore = calculateTotalEnergy();
    
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
    
    // Calculate and store time step metrics
    double energyAfter = calculateTotalEnergy();
    double energyError = std::abs((energyAfter - energyBefore) / energyBefore);
    timeStepHistory.push_back({timeStepSize, energyError, t});
    
    // Update simulation time
    t += timeStepSize;
    
    // Clean up
    delete[] force0;
    delete[] force1;
    delete[] force2;
}

/**
 * Check if simulation has been completed.
 */
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

// void NBodySimulation::printSummary () {
//   std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
//   std::cout << "Position of first remaining object: "
//             << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
// }

void NBodySimulation::printSummary() {
  std::cout << "Simulation Summary:" << std::endl;
  std::cout << "==================" << std::endl;
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
            
  if (!timeStepHistory.empty()) {
      auto bestMetrics = std::min_element(
          timeStepHistory.begin(),
          timeStepHistory.end(),
          [](const TimeStepMetrics& a, const TimeStepMetrics& b) -> bool {
              return a.energyError < b.energyError;
          }
      );
      
      // Calculate average energy error
      double avgError = 0.0;
      for (const auto& metrics : timeStepHistory) {
          avgError += metrics.energyError;
      }
      avgError /= timeStepHistory.size();
      
      std::cout << "\nTime Step Analysis:" << std::endl;
      std::cout << "-------------------" << std::endl;
      std::cout << "Initial time step: " << timeStepHistory.front().timeStep << std::endl;
      std::cout << "Best time step found: " << bestMetrics->timeStep << std::endl;
      std::cout << "Best energy error: " << bestMetrics->energyError << std::endl;
      std::cout << "Average energy error: " << avgError << std::endl;
      std::cout << "Time of best performance: " << bestMetrics->time << std::endl;
      
      // Provide recommendations based on energy error
      std::cout << "\nRecommendation:" << std::endl;
      if (bestMetrics->energyError < 1e-6) {
          std::cout << "Time step performance is excellent - current time step is well-suited for this simulation" << std::endl;
      } else if (bestMetrics->energyError < 1e-4) {
          std::cout << "Time step performance is good - current time step provides good accuracy" << std::endl;
      } else {
          std::cout << "Consider using a smaller time step (try " 
                    << bestMetrics->timeStep * 0.5 
                    << ") for better accuracy" << std::endl;
      }
  }
}