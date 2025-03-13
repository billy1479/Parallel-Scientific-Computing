#include "NBodySimulation.h"
#include <iomanip>

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {};

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
              << "+ One body moving away from the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.1  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body orbiting another" << std::endl
              << "    0.1  100.0  0.001    0.0 0.0 0.0  0.0 0.0 0.0 1.0     1.0 0.0 0.0  0.0 1.0 0.0 0.1" << std::endl
              << "+ Two bodies spiralling" << std::endl
              << "    0.1  100.0  0.001    0.0 0.0 0.0  0.0 0.5 0.0 1.0     1.0 0.0 0.0  0.0 -0.5 0.0 1.0" << std::endl
              << "+ Three-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup (Lagrange point approximation)" << std::endl
              << "    0.1  10.0  0.000001   -0.5 -0.288675 0.0  0.0 0.759836 0.0  1.0   0.5 -0.288675 0.0  -0.658037 -0.379918 0.0  1.0   0.0 0.577350 0.0  0.658037 -0.379918 0.0  1.0 " << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.1  10.0  0.00001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
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

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
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

  return (x[j][direction]-x[i][direction]) * mass[i]*mass[j] / distance3;
}

void NBodySimulation::updateBody () {
  
  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies]();
  double* force1 = new double[NumberOfBodies]();
  double* force2 = new double[NumberOfBodies]();

  if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

  for (int i=0; i<NumberOfBodies; i++) {
          force0[i] = 0.0;
          force1[i] = 0.0;
          force2[i] = 0.0;
  }   

  for (int i=0; i<NumberOfBodies; i++) {
	  for (int j=i+1; j<NumberOfBodies; j++) {
		  if(i!=j){
			  // x,y,z forces acting on particle i.
			  force0[i] += force_calculation(i,j,0);
			  force1[i] += force_calculation(i,j,1);
			  force2[i] += force_calculation(i,j,2);
			   // x,y,z symmetric forces acting on particle j.
			  force0[j] -= force_calculation(i,j,0);
			  force1[j] -= force_calculation(i,j,1);
			  force2[j] -= force_calculation(i,j,2);

		  }
	  }
  }
  for (int i=0; i < NumberOfBodies; i++){
	  x[i][0] = x[i][0] + timeStepSize * v[i][0];
	  x[i][1] = x[i][1] + timeStepSize * v[i][1];
	  x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
  for (int i=0; i < NumberOfBodies; i++){
	  v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
	  v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
	  v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

	  maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
   }
   t += timeStepSize;

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

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}

// Total energy calculation method

double NBodySimulation::calculateTotalEnergy(bool mode) {
  // Calculate total energy as the sum of kinetic and potential energy
  double totalEnergy = 0.0;
  double kineticEnergy = 0.0;
  double potentialEnergy = 0.0;
  
  // Calculate kinetic energy: sum of 0.5 * mass * velocity^2 for each body
  for (int i = 0; i < NumberOfBodies; i++) {
    // Calculate velocity magnitude squared (v^2 = vx^2 + vy^2 + vz^2)
    double velocitySquared = 
      v[i][0] * v[i][0] + 
      v[i][1] * v[i][1] + 
      v[i][2] * v[i][2];
    
    // Add kinetic energy contribution from this body
    kineticEnergy += 0.5 * mass[i] * velocitySquared;
  }
  
  // Calculate potential energy: sum of -G * m_i * m_j / r_ij for each unique pair
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int j = i + 1; j < NumberOfBodies; j++) {
      // Calculate Euclidean distance between bodies i and j
      double distance = sqrt(
        (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
        (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
        (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])
      );
      
      // Add potential energy contribution from this pair
      // The negative sign is because gravitational potential energy is negative
      potentialEnergy -= (mass[i] * mass[j]) / distance;
    }
  }
  
  // Total energy is the sum of kinetic and potential energy
  totalEnergy = kineticEnergy + potentialEnergy;

  if (mode) {
    std::cout << "Initial energy: " << totalEnergy << std::endl;
  } else {
    std::cout << "Final energy: " << totalEnergy << std::endl;
  }

  return totalEnergy;
}

// Outputs stats to csv file for graph

double NBodySimulation::findLargestStableTimeStep(double initialTimeStep, double maxTimeStep, double maxTime, double energyThreshold) {
  std::cout << "Finding largest stable time step..." << std::endl;
  
  double stableTimeStep = initialTimeStep;
  double testTimeStep = initialTimeStep;
  
  // Save original time step
  double originalTimeStep = timeStepSize;
  double originalTime = t;
  
  // Calculate initial energy
  double initialEnergy = calculateTotalEnergy(true);
  
  // Vectors to store backup state
  std::vector<std::vector<double>> xBackup(NumberOfBodies, std::vector<double>(3));
  std::vector<std::vector<double>> vBackup(NumberOfBodies, std::vector<double>(3));
  
  // Create CSV file for energy drift data
  std::ofstream energyFile("energy_stability.csv");
  energyFile << "TimeStep,SimulationTime,Energy,EnergyChange(%),MinDistance,MaxVelocity,Stable" << std::endl;
  
  // Binary search approach for stability threshold
  double lowerBound = initialTimeStep;
  double upperBound = maxTimeStep;
  
  while ((upperBound - lowerBound) / lowerBound > 0.01) { // Continue until 1% precision
      testTimeStep = (lowerBound + upperBound) / 2.0;
      std::cout << "Testing time step: " << testTimeStep << std::endl;
      
      // Backup current state
      backupState(xBackup, vBackup);
      
      // Reset time and set test time step
      t = originalTime;
      timeStepSize = testTimeStep;
      
      // Run simulation for a fixed physical time
      bool unstable = false;
      double testTime = 0.0;
      
      while (testTime < maxTime && !unstable) {
          updateBody();
          testTime += timeStepSize;
          
          // Check for numeric instabilities (NaN or Inf)
          for (int i = 0; i < NumberOfBodies; i++) {
              for (int d = 0; d < 3; d++) {
                  if (std::isnan(x[i][d]) || std::isinf(x[i][d]) ||
                      std::isnan(v[i][d]) || std::isinf(v[i][d])) {
                      unstable = true;
                      break;
                  }
              }
              if (unstable) break;
          }
          
          // Skip additional checks if already unstable
          if (unstable) break;
          
          // Record energy data at regular intervals
          if (fmod(testTime, maxTime/10.0) < timeStepSize) {
              double currentEnergy = calculateTotalEnergy(false);
              double energyChange = std::abs((currentEnergy - initialEnergy) / initialEnergy);
              
              // Write data to CSV
              energyFile << testTimeStep << ","
                        << testTime << ","
                        << currentEnergy << ","
                        << energyChange*100 << ","
                        << minDx << ","
                        << maxV << ","
                        << "1" << std::endl;
              
              // Check stability criterion
              if (energyChange > energyThreshold) {
                  unstable = true;
                  std::cout << "Energy drift detected: " << energyChange*100 << "% change" << std::endl;
                  // Mark as unstable in CSV
                  energyFile << testTimeStep << ","
                            << testTime << ","
                            << currentEnergy << ","
                            << energyChange*100 << ","
                            << minDx << ","
                            << maxV << ","
                            << "0" << std::endl;
                  break;
              }
          }
      }
      
      // Update bounds based on stability
      if (unstable) {
          upperBound = testTimeStep;
          std::cout << "Unstable at dt = " << testTimeStep << std::endl;
      } else {
          lowerBound = testTimeStep;
          stableTimeStep = testTimeStep;
          std::cout << "Stable at dt = " << testTimeStep << std::endl;
      }
      
      // Restore original state
      restoreState(xBackup, vBackup);
  }
  
  // Restore original time step
  timeStepSize = originalTimeStep;
  t = originalTime;
  
  // Close CSV file
  energyFile.close();
  std::cout << "Largest stable time step found: " << stableTimeStep << std::endl;
  std::cout << "Energy stability data written to energy_stability.csv" << std::endl;
  return stableTimeStep;
}

void NBodySimulation::backupState(std::vector<std::vector<double>>& xBackup, std::vector<std::vector<double>>& vBackup) {
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          xBackup[i][d] = x[i][d];
          vBackup[i][d] = v[i][d];
      }
  }
}

void NBodySimulation::restoreState(const std::vector<std::vector<double>>& xBackup, const std::vector<std::vector<double>>& vBackup) {
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          x[i][d] = xBackup[i][d];
          v[i][d] = vBackup[i][d];
      }
  }
}

void NBodySimulation::analyzeConvergenceOrder(double referenceTimeStep, double maxTime, 
  int numRefinements, double energyThreshold) {
std::cout << "Analyzing convergence order of the numerical scheme..." << std::endl;

// Save original simulation state
double originalTimeStep = timeStepSize;
double originalTime = t;

// Store initial state to restore between tests
std::vector<std::vector<double>> initialX(NumberOfBodies, std::vector<double>(3));
std::vector<std::vector<double>> initialV(NumberOfBodies, std::vector<double>(3));

// Save initial state
for (int i = 0; i < NumberOfBodies; i++) {
for (int d = 0; d < 3; d++) {
initialX[i][d] = x[i][d];
initialV[i][d] = v[i][d];
}
}

// Calculate "exact" solution first (using very small time step)
double exactTimeStep = referenceTimeStep / 10.0; // Use a more stable reference time step
std::cout << "Computing reference solution with dt = " << exactTimeStep << std::endl;

// Backup energy before reference run
double initialEnergy = calculateTotalEnergy(false);

// Compute reference solution
double exactFinalEnergy = computeNumericalSolution(exactTimeStep, maxTime, energyThreshold);

// Check if reference solution is valid
if (std::isnan(exactFinalEnergy)) {
std::cerr << "Error: Reference solution resulted in NaN. Try a larger reference time step." << std::endl;
exactFinalEnergy = initialEnergy; // Fallback to initial energy as reference
}

std::cout << "Reference energy value: " << exactFinalEnergy << std::endl;

// Create vectors to store results
std::vector<double> timeSteps;
std::vector<double> energyErrors;
std::vector<double> energyValues;
std::vector<double> minDistances;
std::vector<double> maxVelocities;

// Create CSV file for convergence data
std::ofstream convergenceFile("convergence_analysis.csv");
convergenceFile << "TimeStep,EnergyValue,EnergyError,RelativeError,MinDistance,MaxVelocity,Order" << std::endl;

// Test different time steps
for (int i = 0; i < numRefinements; i++) {
// Calculate time step for this refinement (starting with reference and doubling each time)
double testTimeStep = referenceTimeStep * std::pow(2.0, i);
std::cout << "Testing with time step: " << testTimeStep << std::endl;

// Restore initial state
for (int j = 0; j < NumberOfBodies; j++) {
for (int d = 0; d < 3; d++) {
x[j][d] = initialX[j][d];
v[j][d] = initialV[j][d];
}
}

// Reset time
t = originalTime;

// Track min distance and max velocity during run
double runMinDx = std::numeric_limits<double>::max();
double runMaxV = 0.0;

// Compute solution with this time step
double finalEnergy = computeNumericalSolution(testTimeStep, maxTime, energyThreshold);

// Calculate energy error compared to reference solution
double energyError = 0.0;
double relativeError = 0.0;

if (std::isnan(finalEnergy)) {
// If this simulation became unstable, use a large error value
energyError = 1.0; // Some large value to indicate error
relativeError = 1.0;
std::cout << "  Warning: This run resulted in NaN. Setting error to 1.0" << std::endl;
} else {
energyError = std::abs(finalEnergy - exactFinalEnergy);
if (std::abs(exactFinalEnergy) > 1e-10) {
relativeError = energyError / std::abs(exactFinalEnergy);
} else {
// Avoid division by zero or very small numbers
relativeError = energyError;
}
}

// Store results
timeSteps.push_back(testTimeStep);
energyErrors.push_back(relativeError);
energyValues.push_back(finalEnergy);
minDistances.push_back(minDx); // Use the last calculated minDx
maxVelocities.push_back(maxV); // Use the last calculated maxV

// Calculate convergence order (if we have at least 2 points)
double order = 0.0;
if (i > 0) {
// Calculate order between this point and the previous one
order = std::log(energyErrors[i-1] / energyErrors[i]) / std::log(timeSteps[i] / timeSteps[i-1]);
}

// Write to CSV
convergenceFile << std::scientific << std::setprecision(10)
<< testTimeStep << ","
<< finalEnergy << ","
<< energyError << ","
<< relativeError << ","
<< minDx << ","
<< maxV << ","
<< order << std::endl;

std::cout << "  Final energy: " << finalEnergy 
<< ", Error: " << energyError
<< ", Relative error: " << relativeError * 100 << "%" << std::endl;

if (i > 0) {
std::cout << "  Estimated order of convergence: " << order << std::endl;
}
}

// Calculate the overall convergence order using all data points
double overallOrder = calculateConvergenceOrder(energyErrors, timeSteps);
std::cout << "Overall estimated order of convergence: " << overallOrder << std::endl;

// Write overall order to file
convergenceFile << std::endl << "Overall order of convergence," << overallOrder << std::endl;

// Create a separate file with plotting data in a format easy for visualization tools
std::ofstream plotDataFile("convergence_plot_data.csv");
plotDataFile << "TimeStep,RelativeError" << std::endl;
for (size_t i = 0; i < timeSteps.size(); i++) {
plotDataFile << timeSteps[i] << "," << energyErrors[i] << std::endl;
}

// Add reference lines for first and second order convergence
plotDataFile << std::endl << "# Reference lines" << std::endl;
plotDataFile << "TimeStep,FirstOrder,SecondOrder" << std::endl;

// Use the first valid error point to anchor reference lines
double firstOrderConstant = 0.01;  // Default value
double secondOrderConstant = 0.01; // Default value

// Find first valid error measurement
for (size_t i = 0; i < timeSteps.size(); i++) {
if (!std::isnan(energyErrors[i]) && energyErrors[i] > 0) {
firstOrderConstant = energyErrors[i] / timeSteps[i];
secondOrderConstant = energyErrors[i] / (timeSteps[i] * timeSteps[i]);
break;
}
}

// Generate points for reference lines
for (size_t i = 0; i < timeSteps.size(); i++) {
double dt = timeSteps[i];
plotDataFile << dt << "," 
<< firstOrderConstant * dt << "," 
<< secondOrderConstant * dt * dt << std::endl;
}

// Close files
convergenceFile.close();
plotDataFile.close();

// Restore original time step
timeStepSize = originalTimeStep;
t = originalTime;

// Restore initial state
for (int i = 0; i < NumberOfBodies; i++) {
for (int d = 0; d < 3; d++) {
x[i][d] = initialX[i][d];
v[i][d] = initialV[i][d];
}
}

std::cout << "Convergence analysis completed. Data written to convergence_analysis.csv and convergence_plot_data.csv" << std::endl;
}

double NBodySimulation::computeNumericalSolution(double testTimeStep, double maxTime, double energyThreshold) {
// Set the time step
timeStepSize = testTimeStep;

// Calculate initial energy
double initialEnergy = calculateTotalEnergy(false);

// Run simulation until maxTime
double testTime = 0.0;
bool unstable = false;

while (testTime < maxTime && !unstable) {
updateBody();
testTime += timeStepSize;

// Check for numeric instabilities
for (int i = 0; i < NumberOfBodies; i++) {
for (int d = 0; d < 3; d++) {
if (std::isnan(x[i][d]) || std::isinf(x[i][d]) ||
std::isnan(v[i][d]) || std::isinf(v[i][d])) {
unstable = true;
break;
}
}
if (unstable) break;
}

// Check energy conservation
if (testTime >= maxTime * 0.9) { // Only check near the end
double currentEnergy = calculateTotalEnergy(false);
double energyChange = std::abs((currentEnergy - initialEnergy) / initialEnergy);

if (energyChange > energyThreshold) {
unstable = true;
std::cout << "  Energy drift exceeded threshold: " << energyChange*100 << "%" << std::endl;
}
}
}

// If simulation became unstable, return a special value
if (unstable) {
std::cout << "  Warning: Simulation became unstable with dt = " << testTimeStep << std::endl;
return std::numeric_limits<double>::quiet_NaN();
}

// Return the final energy
return calculateTotalEnergy(false);
}

double NBodySimulation::calculateConvergenceOrder(const std::vector<double>& errors, const std::vector<double>& timeSteps) {
// Use least squares fit for log(error) = log(C) + p*log(dt)
// where p is the order of convergence

int n = errors.size();
if (n < 2) return 0.0;

// Take logs of time steps and errors
std::vector<double> logTimeSteps;
std::vector<double> logErrors;

for (int i = 0; i < n; i++) {
// Skip any NaN values (from unstable simulations)
if (std::isnan(errors[i])) continue;

logTimeSteps.push_back(std::log(timeSteps[i]));
logErrors.push_back(std::log(errors[i]));
}

// If we lost too many points due to instability, return 0
if (logTimeSteps.size() < 2) return 0.0;

// Calculate means
double sumLogDt = 0.0, sumLogError = 0.0;
for (size_t i = 0; i < logTimeSteps.size(); i++) {
sumLogDt += logTimeSteps[i];
sumLogError += logErrors[i];
}
double meanLogDt = sumLogDt / logTimeSteps.size();
double meanLogError = sumLogError / logErrors.size();

// Calculate the slope (which is the order of convergence)
double numerator = 0.0, denominator = 0.0;
for (size_t i = 0; i < logTimeSteps.size(); i++) {
numerator += (logTimeSteps[i] - meanLogDt) * (logErrors[i] - meanLogError);
denominator += (logTimeSteps[i] - meanLogDt) * (logTimeSteps[i] - meanLogDt);
}

if (std::abs(denominator) < 1e-10) return 0.0; // Avoid division by zero

// The order is the slope of the line
return numerator / denominator;
}
