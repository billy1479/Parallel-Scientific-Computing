#include "NBodySimulation.h"

NBodySimulation::NBodySimulation () :
t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
x(nullptr), v(nullptr), mass(nullptr),
x_analytical(nullptr), v_analytical(nullptr), useAnalytical(false),
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
  if (x_analytical != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] x_analytical[i];
    delete [] x_analytical;
  }
  if (v_analytical != nullptr) {
    for (int i=0; i<NumberOfBodies; i++)
      delete [] v_analytical[i];
    delete [] v_analytical;
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

void NBodySimulation::setUp(int argc, char** argv, double dt_ref) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  // timeStepSize = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = dt_ref;
  readArgument++;

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

  useAnalytical = (NumberOfBodies <= 10); // Only use analytical for small systems
  if (useAnalytical) {
    initAnalyticalSolution();
    std::cout << "Using analytical solution for error calculation" << std::endl;
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
  
  if (useAnalytical) {
    updateAnalyticalSolution(t);
  }
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


/// Convergence

void NBodySimulation::initAnalyticalSolution() {
  // Only initialize if we're using analytical solution
  if (!useAnalytical) return;
  
  // Allocate memory for analytical solution if not already done
  if (x_analytical == nullptr) {
    x_analytical = new double*[NumberOfBodies];
    v_analytical = new double*[NumberOfBodies];
    for (int i = 0; i < NumberOfBodies; i++) {
      x_analytical[i] = new double[3];
      v_analytical[i] = new double[3];
      
      // Store initial values
      for (int d = 0; d < 3; d++) {
        x_analytical[i][d] = x[i][d];
        v_analytical[i][d] = v[i][d];
      }
    }
  }
}

void NBodySimulation::updateAnalyticalSolution(double time) {
  if (!useAnalytical || x_analytical == nullptr) return;
  
  // For simple test case with one or two particles with constant velocity
  // We can use x(t) = x₀ + v₀*t for each particle independently
  if (NumberOfBodies <= 2) {
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
        // Use initial position and velocity
        x_analytical[i][d] = x_analytical[i][d] + v_analytical[i][d] * timeStepSize;
      }
    }
  }
  
  // For grid setup with particles initially at rest (no-noise scenario)
  // We can still use constant velocity if no forces are acting
  else if (maxV < 1e-10) { // Particles initially at rest
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int d = 0; d < 3; d++) {
        // Since initial velocity is zero, position doesn't change
        // But update anyway based on current analytical velocity
        x_analytical[i][d] = x_analytical[i][d] + v_analytical[i][d] * timeStepSize;
      }
    }
  }
  
  // If we have gravitational forces, we need to update velocities too
  // This is a simple Euler step for the analytical solution using the exact same
  // force calculation as the main simulation
  else {
    // Calculate forces for the analytical solution
    double* force0 = new double[NumberOfBodies]();
    double* force1 = new double[NumberOfBodies]();
    double* force2 = new double[NumberOfBodies]();
    
    for (int i = 0; i < NumberOfBodies; i++) {
      force0[i] = 0.0;
      force1[i] = 0.0;
      force2[i] = 0.0;
    }
    
    // We use a very small time step for the analytical solution
    double analytical_dt = timeStepSize / 100.0;
    
    // Do 100 small steps to get a more accurate analytical solution
    for (int step = 0; step < 100; step++) {
      // Calculate forces between particles
      for (int i = 0; i < NumberOfBodies; i++) {
        for (int j = i+1; j < NumberOfBodies; j++) {
          if (i != j) {
            // Calculate distance between particles
            double distance = sqrt(
              (x_analytical[j][0]-x_analytical[i][0]) * (x_analytical[j][0]-x_analytical[i][0]) +
              (x_analytical[j][1]-x_analytical[i][1]) * (x_analytical[j][1]-x_analytical[i][1]) +
              (x_analytical[j][2]-x_analytical[i][2]) * (x_analytical[j][2]-x_analytical[i][2])
            );
            double distance3 = distance * distance * distance;
            
            if (distance3 < 1e-10) continue; // Avoid division by zero
            
            // Calculate force components
            force0[i] += (x_analytical[j][0]-x_analytical[i][0]) * mass[i]*mass[j] / distance3;
            force1[i] += (x_analytical[j][1]-x_analytical[i][1]) * mass[i]*mass[j] / distance3;
            force2[i] += (x_analytical[j][2]-x_analytical[i][2]) * mass[i]*mass[j] / distance3;
            
            force0[j] -= (x_analytical[j][0]-x_analytical[i][0]) * mass[i]*mass[j] / distance3;
            force1[j] -= (x_analytical[j][1]-x_analytical[i][1]) * mass[i]*mass[j] / distance3;
            force2[j] -= (x_analytical[j][2]-x_analytical[i][2]) * mass[i]*mass[j] / distance3;
          }
        }
      }
      
      // Update positions using current velocities
      for (int i = 0; i < NumberOfBodies; i++) {
        x_analytical[i][0] += v_analytical[i][0] * analytical_dt;
        x_analytical[i][1] += v_analytical[i][1] * analytical_dt;
        x_analytical[i][2] += v_analytical[i][2] * analytical_dt;
      }
      
      // Update velocities using calculated forces
      for (int i = 0; i < NumberOfBodies; i++) {
        v_analytical[i][0] += force0[i] / mass[i] * analytical_dt;
        v_analytical[i][1] += force1[i] / mass[i] * analytical_dt;
        v_analytical[i][2] += force2[i] / mass[i] * analytical_dt;
      }
      
      // Reset forces for next step
      for (int i = 0; i < NumberOfBodies; i++) {
        force0[i] = 0.0;
        force1[i] = 0.0;
        force2[i] = 0.0;
      }
    }
    
    delete[] force0;
    delete[] force1;
    delete[] force2;
  }
}

double NBodySimulation::calculateError() {
  if (!useAnalytical) {
    return -1.0; // Error calculation not available
  }
  
  if (x_analytical == nullptr) {
    initAnalyticalSolution();
    return 0.0; // No error on first step
  }
  
  double totalError = 0.0;
  for (int i = 0; i < NumberOfBodies; i++) {
    double bodyError = 0.0;
    for (int d = 0; d < 3; d++) {
      double diff = x[i][d] - x_analytical[i][d];
      bodyError += diff * diff;
    }
    totalError += sqrt(bodyError);
  }
  
  // Return average error across all bodies
  return totalError / NumberOfBodies;
}

// Function to compute the error between current and reference solution
double NBodySimulation::computeError(double** x_num, double** x_ref, int numBodies) {
  double error = 0.0;

  for (int i = 0; i < numBodies; i++) {
      double dx = x_num[i][0] - x_ref[i][0];
      double dy = x_num[i][1] - x_ref[i][1];
      double dz = x_num[i][2] - x_ref[i][2];

      double distanceError = sqrt(dx * dx + dy * dy + dz * dz);
      error = std::max(error, distanceError);  // Use max norm
  }

  return error;
}
