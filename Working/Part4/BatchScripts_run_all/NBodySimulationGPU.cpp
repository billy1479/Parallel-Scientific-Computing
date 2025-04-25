#include "NBodySimulationGPU.h"
#include "NBodyInputReader.h"
#include <algorithm>
#include <omp.h>

NBodySimulationGPU::NBodySimulationGPU() :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  force0(nullptr), force1(nullptr), force2(nullptr),
  x_flat(nullptr), v_flat(nullptr), force_flat(nullptr),
  timeStepSize(0), maxV(0), minDx(0),
  snapshotCounter(0), timeStepCounter(0) {
}

NBodySimulationGPU::~NBodySimulationGPU() {
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

  freeForceArrays();
  freeGPUArrays();
}

void NBodySimulationGPU::allocateForceArrays() {
  force0 = new double[NumberOfBodies]();
  force1 = new double[NumberOfBodies]();
  force2 = new double[NumberOfBodies]();
}

void NBodySimulationGPU::freeForceArrays() {
  if (force0 != nullptr) delete[] force0;
  if (force1 != nullptr) delete[] force1;
  if (force2 != nullptr) delete[] force2;
}

void NBodySimulationGPU::setupGPUArrays() {
  // Allocate flat arrays for GPU computation
  x_flat = new double[NumberOfBodies * 3]();
  v_flat = new double[NumberOfBodies * 3]();
  force_flat = new double[NumberOfBodies * 3]();

  // Initialize with data from CPU arrays
  convertToFlatArrays();
}

void NBodySimulationGPU::freeGPUArrays() {
  if (x_flat != nullptr) delete[] x_flat;
  if (v_flat != nullptr) delete[] v_flat;
  if (force_flat != nullptr) delete[] force_flat;
}

void NBodySimulationGPU::convertToFlatArrays() {
  // Fill flat arrays with data from 2D arrays
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int d = 0; d < 3; d++) {
      x_flat[i*3 + d] = x[i][d];
      v_flat[i*3 + d] = v[i][d];
    }
  }
}

void NBodySimulationGPU::convertFromFlatArrays() {
  // Update 2D arrays from flat arrays
  for (int i = 0; i < NumberOfBodies; i++) {
    for (int d = 0; d < 3; d++) {
      x[i][d] = x_flat[i*3 + d];
      v[i][d] = v_flat[i*3 + d];
    }
  }

  // Update force arrays
  for (int i = 0; i < NumberOfBodies; i++) {
    force0[i] = force_flat[i*3];
    force1[i] = force_flat[i*3 + 1];
    force2[i] = force_flat[i*3 + 2];
  }
}

void NBodySimulationGPU::checkInput(int argc, char** argv) {
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

void NBodySimulationGPU::setupFromInputFile(const std::string& filename) {
  try {
    SimulationParams params = NBodyInputReader::readFromFile(filename);

    // Set simulation parameters
    tPlotDelta = params.tPlotDelta;
    tFinal = params.tFinal;
    timeStepSize = params.timeStepSize;
    NumberOfBodies = params.numberOfBodies;

    // Allocate memory for 2D arrays
    x = new double*[NumberOfBodies];
    v = new double*[NumberOfBodies];
    mass = new double[NumberOfBodies];

    // Copy data from params to simulation arrays
    for (int i = 0; i < NumberOfBodies; i++) {
      x[i] = new double[3];
      v[i] = new double[3];

      // Positions
      x[i][0] = params.positions[i*3];
      x[i][1] = params.positions[i*3+1];
      x[i][2] = params.positions[i*3+2];

      // Velocities
      v[i][0] = params.velocities[i*3];
      v[i][1] = params.velocities[i*3+1];
      v[i][2] = params.velocities[i*3+2];

      // Mass
      mass[i] = params.masses[i];
    }

    std::cout << "Read setup with " << NumberOfBodies << " bodies from " << filename << std::endl;

    if (tPlotDelta <= 0.0) {
      std::cout << "plotting switched off" << std::endl;
      tPlot = tFinal + 1.0;
    } else {
      std::cout << "plot initial setup plus every " << tPlotDelta
                << " time units" << std::endl;
      tPlot = 0.0;
    }

  } catch (const std::exception& e) {
    std::cerr << "Error reading input file: " << e.what() << std::endl;
    throw;
  }
}

void NBodySimulationGPU::setupFromCommandLine(int argc, char** argv) {
  // Skip checking input if we're using --input option
  if (!(argc >= 3 && std::string(argv[1]) == "--input")) {
    checkInput(argc, argv);
  }

  NumberOfBodies = (argc-4) / 7;

  x = new double*[NumberOfBodies];
  v = new double*[NumberOfBodies];
  mass = new double[NumberOfBodies];

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
      throw std::runtime_error("Invalid mass value");
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

void NBodySimulationGPU::setUp(int argc, char** argv) {
  // Check if using new input file format
  if (argc >= 3 && std::string(argv[1]) == "--input") {
    std::cout << "Setting up from input file: " << argv[2] << std::endl;
    setupFromInputFile(argv[2]);
  } else {
    // Use old command line format
    std::cout << "Setting up from command line arguments" << std::endl;
    setupFromCommandLine(argc, argv);
  }

  // Initialize force arrays
  allocateForceArrays();

  // Initialize GPU arrays
  setupGPUArrays();

  std::cout << "GPU arrays initialized with " << NumberOfBodies << " bodies" << std::endl;
}

double NBodySimulationGPU::force_calculation(int i, int j, int direction) {
  // Euclidean distance
  const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );
  const double distance3 = distance * distance * distance;
  minDx = std::min(minDx, distance);

  return (x[j][direction]-x[i][direction]) * mass[i]*mass[j] / distance3;
}

void NBodySimulationGPU::runSimulation() {
  timeStepCounter = 0;
  maxV = 0.0;
  minDx = std::numeric_limits<double>::max();

  // Begin the GPU computation - create a target data region to persist data on the GPU
  #pragma omp target data map(tofrom: this[:1], x_flat[0:NumberOfBodies*3], v_flat[0:NumberOfBodies*3], \
                               force_flat[0:NumberOfBodies*3]) \
                        map(to: mass[0:NumberOfBodies], timeStepSize)
  {
    while (!hasReachedEnd()) {
      timeStepCounter++;
      
      // Variables for reduction operations
      double local_minDx = std::numeric_limits<double>::max();
      double local_maxV = 0.0;

      // Reset forces
      #pragma omp target teams distribute parallel for
      for (int i = 0; i < NumberOfBodies * 3; i++) {
        force_flat[i] = 0.0;
      }
      
      #pragma omp target teams distribute parallel for reduction(min:local_minDx)
      for (int i = 0; i < NumberOfBodies; i++) {
        // Thread-private arrays for local force accumulation
        double local_force_x = 0.0;
        double local_force_y = 0.0;
        double local_force_z = 0.0;

        for (int j = 0; j < NumberOfBodies; j++) {
          if (i != j) {
            double dx = x_flat[j * 3] - x_flat[i * 3];
            double dy = x_flat[j * 3 + 1] - x_flat[i * 3 + 1];
            double dz = x_flat[j * 3 + 2] - x_flat[i * 3 + 2];

            double distanceSqr = dx * dx + dy * dy + dz * dz;
            double distance = sqrt(distanceSqr);
            double distance3 = distance * distance * distance;

            // Update minimum distance with reduction
            local_minDx = std::min(local_minDx, distance);

            double factor = mass[i] * mass[j] / distance3;

            // Accumulate forces locally
            local_force_x += dx * factor;
            local_force_y += dy * factor;
            local_force_z += dz * factor;
          }
        }

        // Write back accumulated forces to global arrays
        #pragma omp atomic
        force_flat[i * 3] += local_force_x;
        #pragma omp atomic
        force_flat[i * 3 + 1] += local_force_y;
        #pragma omp atomic
        force_flat[i * 3 + 2] += local_force_z;
      }

      // Update positions and velocities
      #pragma omp target teams distribute parallel for reduction(max:local_maxV)
      for (int i = 0; i < NumberOfBodies; i++) {
        x_flat[i * 3] += timeStepSize * v_flat[i * 3];
        x_flat[i * 3 + 1] += timeStepSize * v_flat[i * 3 + 1];
        x_flat[i * 3 + 2] += timeStepSize * v_flat[i * 3 + 2];

        v_flat[i * 3] += timeStepSize * force_flat[i * 3] / mass[i];
        v_flat[i * 3 + 1] += timeStepSize * force_flat[i * 3 + 1] / mass[i];
        v_flat[i * 3 + 2] += timeStepSize * force_flat[i * 3 + 2] / mass[i];

        double v_squared = v_flat[i * 3] * v_flat[i * 3] +
                           v_flat[i * 3 + 1] * v_flat[i * 3 + 1] +
                           v_flat[i * 3 + 2] * v_flat[i * 3 + 2];

        local_maxV = std::max(local_maxV, sqrt(v_squared));
      }

      // Update CPU metrics with the GPU results
      #pragma omp target update from(local_minDx, local_maxV)
      minDx = local_minDx;
      maxV = local_maxV;

      t += timeStepSize;

      // Take a snapshot if needed
      takeSnapshot();
    }
  }

  // Copy updated positions and velocities back to CPU data structures
  convertFromFlatArrays();
}

void NBodySimulationGPU::updateBody() {
  timeStepCounter++;
  maxV = 0.0;
  minDx = std::numeric_limits<double>::max();

  // Handle special case for a single body
  if (NumberOfBodies == 1) {
    minDx = 0;  // No distances to calculate

    // Update positions
    for (int d = 0; d < 3; d++) {
      x[0][d] += timeStepSize * v[0][d];
    }

    // Calculate maxV
    maxV = std::sqrt(v[0][0]*v[0][0] + v[0][1]*v[0][1] + v[0][2]*v[0][2]);

    t += timeStepSize;
    return;
  }

  // Variables for reduction operations
  double local_minDx = std::numeric_limits<double>::max();
  double local_maxV = 0.0;

  // Begin the GPU computation - create a target data region to minimize data transfers
  #pragma omp target data map(tofrom: x_flat[0:NumberOfBodies*3], v_flat[0:NumberOfBodies*3], \
                               force_flat[0:NumberOfBodies*3], local_minDx, local_maxV) \
                        map(to: mass[0:NumberOfBodies], timeStepSize)
  {
    // Reset forces
    #pragma omp target teams distribute parallel for
    for (int i = 0; i < NumberOfBodies * 3; i++) {
      force_flat[i] = 0.0;
    }

    // Compute forces between all pairs of bodies
    // #pragma omp target teams distribute parallel for collapse(2) reduction(min:local_minDx)
    #pragma omp target teams distribute parallel for reduction(min:local_minDx)
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {
        // Calculate distance between bodies i and j
        double dx = x_flat[j*3] - x_flat[i*3];
        double dy = x_flat[j*3 + 1] - x_flat[i*3 + 1];
        double dz = x_flat[j*3 + 2] - x_flat[i*3 + 2];

        double distanceSqr = dx*dx + dy*dy + dz*dz;
        double distance = sqrt(distanceSqr);
        double distance3 = distance * distance * distance;

        // Update minimum distance
        local_minDx = std::min(local_minDx, distance);

        // Calculate force magnitude
        double factor = mass[i] * mass[j] / distance3;

        // Calculate force components
        double fx = dx * factor;
        double fy = dy * factor;
        double fz = dz * factor;

        // Update forces for body i (needs atomic operations to avoid race conditions)
        #pragma omp atomic
        force_flat[i*3] += fx;
        #pragma omp atomic
        force_flat[i*3 + 1] += fy;
        #pragma omp atomic
        force_flat[i*3 + 2] += fz;

        // Update forces for body j (symmetric but opposite)
        #pragma omp atomic
        force_flat[j*3] -= fx;
        #pragma omp atomic
        force_flat[j*3 + 1] -= fy;
        #pragma omp atomic
        force_flat[j*3 + 2] -= fz;
      }
    }

    // Update positions and velocities
    #pragma omp target teams distribute parallel for reduction(max:local_maxV)
    for (int i = 0; i < NumberOfBodies; i++) {
      // Update positions
      x_flat[i*3]     += timeStepSize * v_flat[i*3];
      x_flat[i*3 + 1] += timeStepSize * v_flat[i*3 + 1];
      x_flat[i*3 + 2] += timeStepSize * v_flat[i*3 + 2];

      // Update velocities
      v_flat[i*3]     += timeStepSize * force_flat[i*3] / mass[i];
      v_flat[i*3 + 1] += timeStepSize * force_flat[i*3 + 1] / mass[i];
      v_flat[i*3 + 2] += timeStepSize * force_flat[i*3 + 2] / mass[i];

      // Calculate velocity magnitude for this body
      double v_squared = v_flat[i*3]*v_flat[i*3] +
                        v_flat[i*3+1]*v_flat[i*3+1] +
                        v_flat[i*3+2]*v_flat[i*3+2];

      // Update maximum velocity
      local_maxV = std::max(local_maxV, sqrt(v_squared));
    }
  } // End of target data region

  // Update CPU metrics with the GPU results
  minDx = local_minDx;
  maxV = local_maxV;

  // Copy updated positions and velocities back to CPU data structures
  convertFromFlatArrays();

  t += timeStepSize;
}

bool NBodySimulationGPU::hasReachedEnd() {
  return t > tFinal;
}

void NBodySimulationGPU::takeSnapshot() {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}

void NBodySimulationGPU::openParaviewVideoFile() {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulationGPU::closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulationGPU::printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out(filename.str().c_str());
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

void NBodySimulationGPU::printSnapshotSummary() {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulationGPU::printSummary() {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}
