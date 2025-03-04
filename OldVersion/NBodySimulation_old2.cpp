#include "NBodySimulation.h"

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

  std::cout << "Original time step: " << timeStepSize << std::endl;

  timeStepSize = findStableTimeStep();

  std::cout << "Largest time step result: " << timeStepSize << std::endl;

  measureConvergenceOrder();
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


// Time step methods

double NBodySimulation::calculateTotalEnergy() {
  double energy = 0.0;
  
  // Kinetic energy: sum of 0.5 * m * v^2 for each body
  for (int i = 0; i < NumberOfBodies; i++) {
      double speedSquared = v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2];
      energy += 0.5 * mass[i] * speedSquared;
  }
  
  // Potential energy: sum of -G * m1 * m2 / r for each pair
  // Note: G = 1 in the simulation units
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i + 1; j < NumberOfBodies; j++) {
          double dx = x[j][0] - x[i][0];
          double dy = x[j][1] - x[i][1];
          double dz = x[j][2] - x[i][2];
          double distance = sqrt(dx*dx + dy*dy + dz*dz);
          
          // Skip if distance is zero (prevents division by zero)
          if (distance > 0) {
              energy -= (mass[i] * mass[j]) / distance;
          }
      }
  }
  
  return energy;
}

double NBodySimulation::findStableTimeStep() {
  std::cout << "\n=== FINDING STABLE TIME STEP ===" << std::endl;
  
  // Step 1: Analyze the physical system
  // -------------------------------------------------------------
  
  // Find minimum distance between any two bodies
  double minDistance = std::numeric_limits<double>::max();
  
  // Find maximum velocity of any body
  double maxVelocity = 0.0;
  
  // Find minimum orbital period between any pair of bodies
  double minOrbitalPeriod = std::numeric_limits<double>::max();
  
  // Calculate force-related quantities
  double maxAcceleration = 0.0;
  
  for (int i = 0; i < NumberOfBodies; i++) {
      // Check velocity of each body
      double velocity = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      maxVelocity = std::max(maxVelocity, velocity);
      
      for (int j = i + 1; j < NumberOfBodies; j++) {
          // Calculate distance between bodies i and j
          double dx = x[j][0] - x[i][0];
          double dy = x[j][1] - x[i][1];
          double dz = x[j][2] - x[i][2];
          double distance = sqrt(dx*dx + dy*dy + dz*dz);
          
          minDistance = std::min(minDistance, distance);
          
          // Only calculate orbital period for non-zero distances
          if (distance > 1e-10) {
              double period = calculateOrbitalPeriod(i, j);
              if (period < minOrbitalPeriod) {
                  minOrbitalPeriod = period;
              }
              
              // Calculate acceleration due to this pair
              double force = mass[i] * mass[j] / (distance * distance);
              double accelI = force / mass[i];
              double accelJ = force / mass[j];
              maxAcceleration = std::max(maxAcceleration, std::max(accelI, accelJ));
          }
      }
  }
  
  // Print physical system properties
  std::cout << "System properties:" << std::endl;
  std::cout << "  * Number of bodies: " << NumberOfBodies << std::endl;
  std::cout << "  * Minimum distance: " << minDistance << std::endl;
  std::cout << "  * Maximum velocity: " << maxVelocity << std::endl;
  std::cout << "  * Minimum orbital period: " << minOrbitalPeriod << std::endl;
  std::cout << "  * Maximum acceleration: " << maxAcceleration << std::endl;
  
  // Step 2: Calculate time step candidates based on different criteria
  // -------------------------------------------------------------
  
  // Initialize with a reasonable default
  double timeStep = tFinal / 1000.0;  // 1/1000 of the total simulation time
  
  // Time step based on distance and velocity (CFL-like condition)
  double dtCFL = std::numeric_limits<double>::max();
  if (maxVelocity > 1e-10) {
      dtCFL = 0.1 * minDistance / maxVelocity;  // Safety factor of 0.1
  }
  
  // Time step based on acceleration
  double dtAccel = std::numeric_limits<double>::max();
  if (maxAcceleration > 1e-10) {
      dtAccel = 0.1 * sqrt(minDistance / maxAcceleration);  // Safety factor of 0.1
  }
  
  // Time step based on orbital period (typically need ~20-100 steps per orbit)
  double dtOrbital = std::numeric_limits<double>::max();
  if (minOrbitalPeriod < std::numeric_limits<double>::max()) {
      dtOrbital = minOrbitalPeriod / 40.0;  // 40 steps per shortest orbit
  }
  
  // Print time step candidates
  std::cout << "\nTime step candidates:" << std::endl;
  std::cout << "  * CFL condition: " << dtCFL << std::endl;
  std::cout << "  * Acceleration-based: " << dtAccel << std::endl;
  std::cout << "  * Orbital period-based: " << dtOrbital << std::endl;
  
  // Step 3: Take the most restrictive (smallest) time step among the candidates
  // -------------------------------------------------------------
  
  timeStep = std::min(timeStep, dtCFL);
  timeStep = std::min(timeStep, dtAccel);
  timeStep = std::min(timeStep, dtOrbital);
  
  // Make sure the time step is not too small
  if (timeStep < 1e-10) {
      timeStep = 1e-10;
      std::cout << "Warning: Calculated time step is extremely small. Using minimum value." << std::endl;
  }
  
  // Step 4: Test this time step for energy conservation
  // -------------------------------------------------------------
  
  // Save the initial state
  double** xSave = new double*[NumberOfBodies];
  double** vSave = new double*[NumberOfBodies];
  
  for (int i = 0; i < NumberOfBodies; i++) {
      xSave[i] = new double[3];
      vSave[i] = new double[3];
      
      for (int d = 0; d < 3; d++) {
          xSave[i][d] = x[i][d];
          vSave[i][d] = v[i][d];
      }
  }
  
  // Calculate initial energy
  double initialEnergy = calculateTotalEnergy();
  std::cout << "\nInitial energy: " << initialEnergy << std::endl;
  
  // Choose number of test steps (~1 orbital period)
  int testSteps = 100;
  if (minOrbitalPeriod < std::numeric_limits<double>::max()) {
      testSteps = static_cast<int>(minOrbitalPeriod / timeStep);
      testSteps = std::min(testSteps, 1000);  // Cap to avoid excessive computation
  }
  
  // Test simulation with the calculated time step
  double tSave = t;  // Save current time
  double origTimeStepSize = timeStepSize;  // Save original time step
  
  timeStepSize = timeStep;  // Set the candidate time step
  
  std::cout << "Testing time step " << timeStep << " for " << testSteps << " steps..." << std::endl;
  
  bool stable = true;
  double maxEnergyError = 0.0;
  
  for (int step = 0; step < testSteps && stable; step++) {
      updateBody();  // Advance the simulation
      
      // Check energy conservation periodically
      if (step % 10 == 0 || step == testSteps - 1) {
          double currentEnergy = calculateTotalEnergy();
          double relativeError = std::abs((currentEnergy - initialEnergy) / initialEnergy);
          
          maxEnergyError = std::max(maxEnergyError, relativeError);
          
          // Consider unstable if energy error exceeds threshold or becomes NaN
          if (relativeError > 0.01 || std::isnan(relativeError)) {
              stable = false;
              std::cout << "Instability detected at step " << step 
                        << ": energy error = " << relativeError * 100 << "%" << std::endl;
          }
          
          // Also check for NaN in positions and velocities
          for (int i = 0; i < NumberOfBodies; i++) {
              for (int d = 0; d < 3; d++) {
                  if (std::isnan(x[i][d]) || std::isnan(v[i][d])) {
                      stable = false;
                      std::cout << "NaN values detected at step " << step << std::endl;
                      break;
                  }
              }
              if (!stable) break;
          }
      }
  }
  
  std::cout << "Max energy error: " << maxEnergyError * 100 << "%" << std::endl;
  
  // Step 5: If not stable, try smaller time steps
  // -------------------------------------------------------------
  
  if (!stable) {
      // Try a range of smaller time steps
      std::vector<double> testSteps = {timeStep * 0.5, timeStep * 0.2, timeStep * 0.1, timeStep * 0.05};
      
      for (double testTimeStep : testSteps) {
          // Reset simulation state
          t = tSave;
          for (int i = 0; i < NumberOfBodies; i++) {
              for (int d = 0; d < 3; d++) {
                  x[i][d] = xSave[i][d];
                  v[i][d] = vSave[i][d];
              }
          }
          
          // Set new test time step
          timeStepSize = testTimeStep;
          
          // Test this smaller time step
          std::cout << "Testing smaller time step " << testTimeStep << "..." << std::endl;
          
          stable = true;
          maxEnergyError = 0.0;
          int numTestSteps = std::min(static_cast<int>(minOrbitalPeriod / testTimeStep), 1000);
          
          for (int step = 0; step < numTestSteps && stable; step++) {
              updateBody();
              
              // Check energy every 10 steps
              if (step % 10 == 0 || step == numTestSteps - 1) {
                  double currentEnergy = calculateTotalEnergy();
                  double relativeError = std::abs((currentEnergy - initialEnergy) / initialEnergy);
                  
                  maxEnergyError = std::max(maxEnergyError, relativeError);
                  
                  if (relativeError > 0.01 || std::isnan(relativeError)) {
                      stable = false;
                  }
                  
                  // Check for NaN values
                  for (int i = 0; i < NumberOfBodies; i++) {
                      for (int d = 0; d < 3; d++) {
                          if (std::isnan(x[i][d]) || std::isnan(v[i][d])) {
                              stable = false;
                              break;
                          }
                      }
                      if (!stable) break;
                  }
              }
          }
          
          std::cout << "Max energy error: " << maxEnergyError * 100 << "%" << std::endl;
          
          if (stable) {
              timeStep = testTimeStep;
              std::cout << "Found stable time step: " << timeStep << std::endl;
              break;
          }
      }
  }
  
  // Step 6: Restore original state and finalize
  // -------------------------------------------------------------
  
  // Restore the original simulation state
  t = tSave;
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          x[i][d] = xSave[i][d];
          v[i][d] = vSave[i][d];
      }
  }
  
  // Clean up temporary arrays
  for (int i = 0; i < NumberOfBodies; i++) {
      delete[] xSave[i];
      delete[] vSave[i];
  }
  delete[] xSave;
  delete[] vSave;
  
  // Apply a safety factor for the final recommendation
  double recommendedTimeStep = timeStep * 0.9;  // 10% safety margin
  
  // Restore original time step (caller will decide whether to use the recommended one)
  timeStepSize = origTimeStepSize;
  
  std::cout << "\n=== FINAL RECOMMENDATION ===" << std::endl;
  std::cout << "Largest stable time step found: " << timeStep << std::endl;
  std::cout << "Recommended time step (with safety margin): " << recommendedTimeStep << std::endl;
  
  if (minOrbitalPeriod < std::numeric_limits<double>::max()) {
      std::cout << "This gives approximately " 
                << (minOrbitalPeriod / recommendedTimeStep) 
                << " steps per shortest orbital period" << std::endl;
  }
  
  return recommendedTimeStep;
}

double NBodySimulation::calculateOrbitalPeriod(int i, int j) {
  // Calculate distance between bodies
  double dx = x[j][0] - x[i][0];
  double dy = x[j][1] - x[i][1];
  double dz = x[j][2] - x[i][2];
  double distance = sqrt(dx*dx + dy*dy + dz*dz);
  
  // Calculate relative velocity magnitude
  double dvx = v[j][0] - v[i][0];
  double dvy = v[j][1] - v[i][1];
  double dvz = v[j][2] - v[i][2];
  double relVelocity = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
  
  // For a two-body system, orbital period is approximately:
  // T = 2π * sqrt(r³/(G*(m1+m2)))
  // Where G = 1 in simulation units
  double combinedMass = mass[i] + mass[j];
  
  // If bodies aren't moving relative to each other, return a very large period
  if (relVelocity < 1e-10 || distance < 1e-10 || combinedMass < 1e-10) {
      return std::numeric_limits<double>::max();
  }
  
  // Kepler's third law for period
  return 2.0 * M_PI * sqrt(pow(distance, 3) / combinedMass);
}

// Convergence order

void NBodySimulation::measureConvergenceOrder() {
  std::cout << "\n=== CONVERGENCE ORDER ANALYSIS ===" << std::endl;
  
  // Step 1: Choose a reference time step that is known to be stable
  // This should be small enough to produce an accurate solution
  double baseTimeStep = 1.0e-6;  // Start with a very small time step
  // double baseTimeStep = 1;
  // Use the stable time step finder if available
  // baseTimeStep = findStableTimeStep();
  
  std::cout << "Using base time step: " << baseTimeStep << std::endl;
  
  // Step 2: Save the initial state of the system
  std::vector<std::vector<double>> initialX(NumberOfBodies, std::vector<double>(3));
  std::vector<std::vector<double>> initialV(NumberOfBodies, std::vector<double>(3));
  
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          initialX[i][d] = x[i][d];
          initialV[i][d] = v[i][d];
      }
  }
  
  // Step 3: Define test time steps as multiples of the base time step
  // Use a series of time steps with factor of 2 between each
  std::vector<double> timeSteps;
  timeSteps.push_back(baseTimeStep);
  timeSteps.push_back(baseTimeStep * 2);
  timeSteps.push_back(baseTimeStep * 4);
  timeSteps.push_back(baseTimeStep * 8);
  timeSteps.push_back(baseTimeStep * 16);
  
  // Step 4: Choose a fixed end time for all simulations
  // This should be long enough to see meaningful dynamics but not too computationally expensive
  double testDuration = baseTimeStep * 100;  // Simulate for 100 steps at the base time step
  
  // Step 5: Generate reference solution with the smallest time step
  std::cout << "Generating reference solution with dt = " << baseTimeStep << std::endl;
  
  // Store the original time step and time
  double originalTimeStep = timeStepSize;
  double originalTime = t;
  
  // Set the smallest time step
  timeStepSize = baseTimeStep;
  t = 0.0;
  
  // Run simulation until reaching the test duration
  int steps = static_cast<int>(testDuration / baseTimeStep);
  for (int i = 0; i < steps; i++) {
      updateBody();
  }
  
  // Save the reference solution
  std::vector<std::vector<double>> referenceX(NumberOfBodies, std::vector<double>(3));
  std::vector<std::vector<double>> referenceV(NumberOfBodies, std::vector<double>(3));
  
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          referenceX[i][d] = x[i][d];
          referenceV[i][d] = v[i][d];
      }
  }
  
  // Step 6: Run simulations with increasingly larger time steps
  std::vector<double> errors;
  
  for (size_t i = 1; i < timeSteps.size(); i++) {
      double currentTimeStep = timeSteps[i];
      
      std::cout << "Testing with dt = " << currentTimeStep << std::endl;
      
      // Reset to initial state
      for (int j = 0; j < NumberOfBodies; j++) {
          for (int d = 0; d < 3; d++) {
              x[j][d] = initialX[j][d];
              v[j][d] = initialV[j][d];
          }
      }
      
      // Set current time step and reset time
      timeStepSize = currentTimeStep;
      t = 0.0;
      
      // Calculate number of steps needed to reach the same end time
      int currentSteps = static_cast<int>(testDuration / currentTimeStep);
      
      // Run simulation
      for (int step = 0; step < currentSteps; step++) {
          updateBody();
      }
      
      // Calculate error compared to reference solution
      double sumSquaredError = 0.0;
      for (int j = 0; j < NumberOfBodies; j++) {
          for (int d = 0; d < 3; d++) {
              double positionError = x[j][d] - referenceX[j][d];
              double velocityError = v[j][d] - referenceV[j][d];
              
              sumSquaredError += positionError * positionError;
              sumSquaredError += velocityError * velocityError;
          }
      }
      
      // RMS error (Root Mean Square)
      double rmsError = std::sqrt(sumSquaredError / (NumberOfBodies * 6));  // 6 = 3 position + 3 velocity components
      errors.push_back(rmsError);
      
      std::cout << "  Time step: " << currentTimeStep 
                << ", Error: " << rmsError << std::endl;
  }
  
  // Step 7: Calculate convergence order using consecutive pairs of time steps
  std::cout << "\nConvergence order estimates:" << std::endl;
  
  for (size_t i = 0; i < errors.size() - 1; i++) {
      double timeStepRatio = timeSteps[i+2] / timeSteps[i+1];  // Ratio between consecutive time steps
      double errorRatio = errors[i+1] / errors[i];  // Ratio between consecutive errors
      
      // Convergence order = log(error ratio) / log(time step ratio)
      double convergenceOrder = std::log(errorRatio) / std::log(timeStepRatio);
      
      std::cout << "  Between dt = " << timeSteps[i+1] << " and dt = " << timeSteps[i+2]
                << ": Order = " << convergenceOrder << std::endl;
  }
  
  // Step 8: Calculate overall convergence order using least squares fit
  double sumLogDt = 0.0;
  double sumLogError = 0.0;
  double sumLogDtSquared = 0.0;
  double sumLogDtLogError = 0.0;
  int n = errors.size();
  
  for (int i = 0; i < n; i++) {
      double logDt = std::log(timeSteps[i+1]);
      double logError = std::log(errors[i]);
      
      sumLogDt += logDt;
      sumLogError += logError;
      sumLogDtSquared += logDt * logDt;
      sumLogDtLogError += logDt * logError;
  }
  
  // Linear regression to find slope (which equals the convergence order)
  double overallOrder = (n * sumLogDtLogError - sumLogDt * sumLogError) / 
                       (n * sumLogDtSquared - sumLogDt * sumLogDt);
  
  std::cout << "\nOverall estimated convergence order: " << overallOrder << std::endl;
  
  // For a symplectic Euler method, we expect an order of approximately 1
  std::cout << "Expected theoretical order for symplectic Euler: 1.0" << std::endl;
  
  // If using a higher-order method, adjust the expected value accordingly
  // For velocity Verlet / leapfrog: 2.0
  // For 4th-order Runge-Kutta: 4.0
  
  // Step 9: Restore original state
  timeStepSize = originalTimeStep;
  t = originalTime;
  
  for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
          x[i][d] = initialX[i][d];
          v[i][d] = initialV[i][d];
      }
  }
  
  std::cout << "Convergence analysis complete, original simulation state restored." << std::endl;
}
