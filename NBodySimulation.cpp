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

  calculateStableTimeStep();
  
  // Calculate initial energy for stability monitoring
  initialEnergy = calculateTotalEnergy();
  std::cout << "Initial total energy: " << initialEnergy << std::endl;
  std::cout << "Calculated stable timestep: " << timeStepSize << std::endl;
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

// void NBodySimulation::updateBody () {

//   timeStepCounter++;
//   maxV   = 0.0;
//   minDx  = std::numeric_limits<double>::max();

//   // force0 = force along x direction
//   // force1 = force along y direction
//   // force2 = force along z direction
//   double* force0 = new double[NumberOfBodies]();
//   double* force1 = new double[NumberOfBodies]();
//   double* force2 = new double[NumberOfBodies]();

//   if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

//   for (int i=0; i<NumberOfBodies; i++) {
//           force0[i] = 0.0;
//           force1[i] = 0.0;
//           force2[i] = 0.0;
//   }   

//   for (int i=0; i<NumberOfBodies; i++) {
// 	  for (int j=i+1; j<NumberOfBodies; j++) {
// 		  if(i!=j){
// 			  // x,y,z forces acting on particle i.
// 			  force0[i] += force_calculation(i,j,0);
// 			  force1[i] += force_calculation(i,j,1);
// 			  force2[i] += force_calculation(i,j,2);
// 			   // x,y,z symmetric forces acting on particle j.
// 			  force0[j] -= force_calculation(i,j,0);
// 			  force1[j] -= force_calculation(i,j,1);
// 			  force2[j] -= force_calculation(i,j,2);

// 		  }
// 	  }
//   }
//   for (int i=0; i < NumberOfBodies; i++){
// 	  x[i][0] = x[i][0] + timeStepSize * v[i][0];
// 	  x[i][1] = x[i][1] + timeStepSize * v[i][1];
// 	  x[i][2] = x[i][2] + timeStepSize * v[i][2];
//   }
//   for (int i=0; i < NumberOfBodies; i++){
// 	  v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
// 	  v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
// 	  v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

// 	  maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
//    }
//    t += timeStepSize;

//   delete[] force0;
//   delete[] force1;
//   delete[] force2;
// }


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
  
  // Update positions and then velocities (explicit Euler)
  for (int i=0; i < NumberOfBodies; i++){
	  x[i][0] = x[i][0] + timeStepSize * v[i][0];
	  x[i][1] = x[i][1] + timeStepSize * v[i][1];
	  x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
  
  for (int i=0; i < NumberOfBodies; i++){
	  v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
	  v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
	  v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];
      
      // Calculate max velocity for stability analysis
	  maxV = std::max(maxV, std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]));
  }
  
  // Recalculate minimum distance between all particles for adaptive timestep
  for (int i=0; i < NumberOfBodies; i++) {
    for (int j=i+1; j < NumberOfBodies; j++) {
      const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
      minDx = std::min(minDx, distance);
    }
  }
  
  // The explicit Euler method is inherently unstable for this type of simulation
  // Let's implement a fixed timestep that's small enough for stability
  const double fixedTimeStep = 1e-6; // Very small fixed timestep
  
  // Calculate current energy to monitor stability
  currentEnergy = calculateTotalEnergy();
  
  // Only track energy error for monitoring, don't use for adaptive timestepping
  double energyError = 0.0;
  if (std::abs(initialEnergy) > 1e-10) {
    energyError = std::abs((currentEnergy - initialEnergy) / initialEnergy);
    
    // Print energy error but don't adjust timestep
    if (energyError > 0.001) {
      std::cout << "Energy error " << energyError * 100 << "% - monitoring" << std::endl;
    }
  }
  
  // Use fixed timestep instead of adaptive
  timeStepSize = fixedTimeStep;
  
  // Update the time step size
  // timeStepSize = suggestedTimeStep;
  
  // Update simulation time
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

// Time stepping

double NBodySimulation::calculateTotalEnergy() {
  double kineticEnergy = 0.0;
  double potentialEnergy = 0.0;
  
  // Calculate kinetic energy
  for (int i = 0; i < NumberOfBodies; i++) {
    double v_squared = v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
    kineticEnergy += 0.5 * mass[i] * v_squared;
  }
  
  // Calculate potential energy
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int j = i+1; j < NumberOfBodies; j++) {
      const double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
      
      // Gravitational potential energy (G=1 in this code)
      potentialEnergy -= (mass[i] * mass[j]) / distance;
    }
  }
  
  return kineticEnergy + potentialEnergy;
}

void NBodySimulation::calculateStableTimeStep() {
  // Calculate maximum velocity
  maxV = 0.0;
  for (int i = 0; i < NumberOfBodies; i++) {
    double v_mag = sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    maxV = std::max(maxV, v_mag);
  }
  
  // Calculate minimum distance between bodies
  minDx = std::numeric_limits<double>::max();
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int j = i+1; j < NumberOfBodies; j++) {
      double distance = sqrt(
        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
      );
      
      if (distance > 0) {
        minDx = std::min(minDx, distance);
      }
    }
  }
  
  // If there's only one body or all bodies are at the same position
  if (NumberOfBodies <= 1 || minDx == std::numeric_limits<double>::max()) {
    minDx = 1.0;
  }
  
  // Calculate maximum acceleration
  double maxAccel = 0.0;
  for (int i = 0; i < NumberOfBodies; i++) {
    double force[3] = {0.0, 0.0, 0.0};
    
    for (int j = 0; j < NumberOfBodies; j++) {
      if (i != j) {
        double dx = x[j][0] - x[i][0];
        double dy = x[j][1] - x[i][1];
        double dz = x[j][2] - x[i][2];
        
        double distance = sqrt(dx*dx + dy*dy + dz*dz);
        
        if (distance > 1e-10) {
          double distance3 = distance * distance * distance;
          
          // Calculate gravitational force (G=1)
          force[0] += dx * mass[j] / distance3;
          force[1] += dy * mass[j] / distance3;
          force[2] += dz * mass[j] / distance3;
        }
      }
    }
    
    // Calculate acceleration magnitude (F=ma, so a=F/m)
    double accel = sqrt(
      force[0]*force[0] + 
      force[1]*force[1] + 
      force[2]*force[2]
    );
    
    maxAccel = std::max(maxAccel, accel);
  }
  
  // Prevent division by zero
  if (maxAccel < 1e-10) maxAccel = 1e-10;
  if (maxV < 1e-10) maxV = 1e-10;
  
  // Calculate stable timestep based on different criteria
  
  // 1. Position criterion: dt < 0.1 * minDx / maxV
  double dt_pos = 0.1 * minDx / maxV;
  
  // 2. Velocity criterion: dt < 0.1 * sqrt(minDx / maxAccel)
  double dt_vel = 0.1 * sqrt(minDx / maxAccel);
  
  // 3. Acceleration criterion: dt < 0.1 * maxV / maxAccel
  double dt_accel = 0.1 * maxV / maxAccel;
  
  // Take the most restrictive timestep
  timeStepSize = std::min(std::min(dt_pos, dt_vel), dt_accel);
  
  // Add a safety factor and ensure minimum/maximum bounds
  timeStepSize *= 0.5; // 50% safety factor
  
  // Apply reasonable bounds
  const double minTimeStep = 1e-10;
  const double maxTimeStep = 0.01;
  timeStepSize = std::max(minTimeStep, std::min(maxTimeStep, timeStepSize));
  
  std::cout << "Stable timestep calculation:" << std::endl;
  std::cout << "  Max velocity: " << maxV << std::endl;
  std::cout << "  Min distance: " << minDx << std::endl;
  std::cout << "  Max acceleration: " << maxAccel << std::endl;
  std::cout << "  Position criterion: " << dt_pos << std::endl;
  std::cout << "  Velocity criterion: " << dt_vel << std::endl;
  std::cout << "  Acceleration criterion: " << dt_accel << std::endl;
}