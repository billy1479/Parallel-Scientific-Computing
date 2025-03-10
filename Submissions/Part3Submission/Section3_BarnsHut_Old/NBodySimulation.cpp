#include "NBodySimulation.h"
#include "OctreeNode.h"
#include <cmath>
#include <limits>

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0), domainSize(0), octree(nullptr) {};

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
  if (octree != nullptr) {
    delete octree;
    octree = nullptr;
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

  domainSize = findDomainSize();

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

// Octree methods
double NBodySimulation::findDomainSize() {
  double maxCoord = 0.0;
  
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int d = 0; d < 3; d++) {
      maxCoord = std::max(maxCoord, std::abs(x[i][d]));
    }
  }
  
  // Add some margin and ensure power of 2 for efficient octree division
  double size = maxCoord * 2.0;
  
  // Round up to power of 2
  return std::pow(2, std::ceil(std::log2(size)));
}

// Octree building method
void NBodySimulation::buildOctree() {
  std::cout << "Starting buildOctree()" << std::endl;
  // Delete old octree if it exists
  if (octree != nullptr) {
    std::cout << "Deleting old octree" << std::endl;
    delete octree;
    octree = nullptr;
  }

  std::cout << "Creating new octree with domainSize = " << domainSize << std::endl;
  
  // Create new octree
  octree = new OctreeNode(0.0, 0.0, 0.0, domainSize / 2.0);

  std::cout << "Inserting " << NumberOfBodies << " bodies into octree" << std::endl;
  
  // Insert all bodies
  for (int i = 0; i < NumberOfBodies; i++) {
    octree->insert(i, x, mass);
  }

  std::cout << "Finished buildOctree()" << std::endl;
}

void NBodySimulation::updateBody() {
  timeStepCounter++;
  maxV = 0.0;
  minDx = std::numeric_limits<double>::max();

  // Build the octree for this time step
  buildOctree();

  // Calculate forces using Barnes-Hut approximation
  double* force0 = new double[NumberOfBodies]();
  double* force1 = new double[NumberOfBodies]();
  double* force2 = new double[NumberOfBodies]();

  // Reset forces
  for (int i = 0; i < NumberOfBodies; i++) {
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;
  }

  // For each body, compute forces using the octree
  for (int i = 0; i < NumberOfBodies; i++) {
    double forces[3] = {0.0, 0.0, 0.0};
    octree->computeForce(i, x, forces, mass);
    
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

  // Update velocities using computed forces
  for (int i = 0; i < NumberOfBodies; i++) {
    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

    // Update maximum velocity
    maxV = std::max(maxV, std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]));
  }

  // Recalculate minimum distance between bodies (optional, for reporting)
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

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
}