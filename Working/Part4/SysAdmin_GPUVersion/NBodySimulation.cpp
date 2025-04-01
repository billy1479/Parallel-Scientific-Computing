#include "NBodySimulation.h"
#include <omp.h>  // OpenMP header
#include <iomanip> // For std::setprecision

// Initialize timing variables
double force_calc_time = 0.0;
double position_update_time = 0.0;
double velocity_update_time = 0.0;
double total_simulation_time = 0.0;
double data_transfer_time = 0.0;  // Track data transfer time to/from GPU

NBodySimulation::NBodySimulation () :
  force_calc_time(0.0), position_update_time(0.0), velocity_update_time(0.0), 
  total_simulation_time(0.0), data_transfer_time(0.0),
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {
  
  // Create and initialize the timing CSV file
  timing_file.open("nbody_timing.csv");
  // Write CSV header
  timing_file << "TimeStep,SimulationTime,TotalElapsedTime,ForceCalcTime,PositionUpdateTime,VelocityUpdateTime,"
              << "DataTransferTime,ForceCalcPercentage,PositionUpdatePercentage,VelocityUpdatePercentage,"
              << "DataTransferPercentage,NumberOfBodies" << std::endl;
};

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
  
  // Print timing summary at the end
  std::cout << "\n=== Timing Summary ===\n";
  std::cout << "Total simulation time: " << total_simulation_time << " seconds\n";
  std::cout << "Force calculation time: " << force_calc_time << " seconds (" 
            << (force_calc_time/total_simulation_time*100) << "%)\n";
  std::cout << "Position update time: " << position_update_time << " seconds ("
            << (position_update_time/total_simulation_time*100) << "%)\n";
  std::cout << "Velocity update time: " << velocity_update_time << " seconds ("
            << (velocity_update_time/total_simulation_time*100) << "%)\n";
  std::cout << "Data transfer time: " << data_transfer_time << " seconds ("
            << (data_transfer_time/total_simulation_time*100) << "%)\n";
  
  // Write final summary to CSV and close file
  timing_file << "SUMMARY,,"
              << total_simulation_time << ","
              << force_calc_time << ","
              << position_update_time << ","
              << velocity_update_time << ","
              << data_transfer_time << ","
              << (force_calc_time/total_simulation_time*100) << ","
              << (position_update_time/total_simulation_time*100) << ","
              << (velocity_update_time/total_simulation_time*100) << ","
              << (data_transfer_time/total_simulation_time*100) << ","
              << NumberOfBodies << std::endl;
  
  // Close the timing file
  timing_file.close();
  
  std::cout << "Timing data has been written to nbody_timing.csv" << std::endl;
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

void NBodySimulation::updateBody() {
  double start_time, end_time, section_start;
  double step_force_time = 0.0;
  double step_position_time = 0.0;
  double step_velocity_time = 0.0;
  double step_data_transfer_time = 0.0;
  double step_total_time = 0.0;
  
  // Start measuring total time for this update
  start_time = omp_get_wtime();

  timeStepCounter++;
  maxV = 0.0;
  minDx = std::numeric_limits<double>::max();

  // Prepare data for GPU processing
  // We need to reorganize our data structure for GPU offloading
  // Instead of an array of pointers, we need flat arrays
  
  // Prepare flat arrays for positions, velocities, and forces
  double* pos_x = new double[NumberOfBodies]; 
  double* pos_y = new double[NumberOfBodies];
  double* pos_z = new double[NumberOfBodies];
  
  double* vel_x = new double[NumberOfBodies];
  double* vel_y = new double[NumberOfBodies];
  double* vel_z = new double[NumberOfBodies];
  
  double* force_x = new double[NumberOfBodies]();  // Initialize to zero
  double* force_y = new double[NumberOfBodies]();
  double* force_z = new double[NumberOfBodies]();
  
  double* mass_flat = new double[NumberOfBodies];
  
  // Start measuring data transfer time (host to device setup)
  section_start = omp_get_wtime();
  
  // Copy data to flat arrays
  for (int i = 0; i < NumberOfBodies; i++) {
    pos_x[i] = x[i][0];
    pos_y[i] = x[i][1];
    pos_z[i] = x[i][2];
    
    vel_x[i] = v[i][0];
    vel_y[i] = v[i][1];
    vel_z[i] = v[i][2];
    
    mass_flat[i] = mass[i];
  }
  
  step_data_transfer_time += omp_get_wtime() - section_start;
  data_transfer_time += step_data_transfer_time;
  
  // Force calculation section - offload to GPU
  section_start = omp_get_wtime();
  
  // Local variable for minimum distance tracking on the device
  double device_min_dx = std::numeric_limits<double>::max();
  
  // Begin target region - offload work to GPU
  #pragma omp target data map(to: pos_x[0:NumberOfBodies], pos_y[0:NumberOfBodies], pos_z[0:NumberOfBodies], \
                              mass_flat[0:NumberOfBodies]) \
                         map(tofrom: force_x[0:NumberOfBodies], force_y[0:NumberOfBodies], force_z[0:NumberOfBodies]) \
                         map(from: device_min_dx)
  {
    // Use a non-collapsed approach for the nested loops
    #pragma omp target teams distribute parallel for reduction(min:device_min_dx)
    for (int i = 0; i < NumberOfBodies; i++) {
      for (int j = i+1; j < NumberOfBodies; j++) {
        // Calculate distance
        const double dx = pos_x[j] - pos_x[i];
        const double dy = pos_y[j] - pos_y[i];
        const double dz = pos_z[j] - pos_z[i];
        
        const double distance = sqrt(dx*dx + dy*dy + dz*dz);
        const double distance3 = distance * distance * distance;
        
        // Update local minimum distance
        device_min_dx = (distance < device_min_dx) ? distance : device_min_dx;
        
        // Calculate force components
        const double force_magnitude = mass_flat[i] * mass_flat[j] / distance3;
        
        // Use atomic operations to avoid race conditions when updating forces
        #pragma omp atomic update
        force_x[i] += dx * force_magnitude;
        
        #pragma omp atomic update
        force_y[i] += dy * force_magnitude;
        
        #pragma omp atomic update
        force_z[i] += dz * force_magnitude;
        
        #pragma omp atomic update
        force_x[j] -= dx * force_magnitude;
        
        #pragma omp atomic update
        force_y[j] -= dy * force_magnitude;
        
        #pragma omp atomic update
        force_z[j] -= dz * force_magnitude;
      }
    }
  } // End target data region
  
  // Update the class member with the device-computed minimum distance
  minDx = std::min(minDx, device_min_dx);
  
  step_force_time = omp_get_wtime() - section_start;
  force_calc_time += step_force_time;
  
  // Position update section
  section_start = omp_get_wtime();
  
  #pragma omp target teams distribute parallel for map(tofrom: pos_x[0:NumberOfBodies], pos_y[0:NumberOfBodies], pos_z[0:NumberOfBodies]) \
                                               map(to: vel_x[0:NumberOfBodies], vel_y[0:NumberOfBodies], vel_z[0:NumberOfBodies])
  for (int i = 0; i < NumberOfBodies; i++) {
    pos_x[i] = pos_x[i] + timeStepSize * vel_x[i];
    pos_y[i] = pos_y[i] + timeStepSize * vel_y[i];
    pos_z[i] = pos_z[i] + timeStepSize * vel_z[i];
  }
  
  step_position_time = omp_get_wtime() - section_start;
  position_update_time += step_position_time;
  
  // Velocity update section
  section_start = omp_get_wtime();
  
  double local_maxV = 0.0;
  
  #pragma omp target teams distribute parallel for map(tofrom: vel_x[0:NumberOfBodies], vel_y[0:NumberOfBodies], vel_z[0:NumberOfBodies]) \
                                               map(to: force_x[0:NumberOfBodies], force_y[0:NumberOfBodies], force_z[0:NumberOfBodies], \
                                                   mass_flat[0:NumberOfBodies]) \
                                               reduction(max:local_maxV)
  for (int i = 0; i < NumberOfBodies; i++) {
    vel_x[i] = vel_x[i] + timeStepSize * force_x[i] / mass_flat[i];
    vel_y[i] = vel_y[i] + timeStepSize * force_y[i] / mass_flat[i];
    vel_z[i] = vel_z[i] + timeStepSize * force_z[i] / mass_flat[i];
    
    double velocity_magnitude = sqrt(vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
    local_maxV = (velocity_magnitude > local_maxV) ? velocity_magnitude : local_maxV;
  }
  
  maxV = local_maxV;
  
  step_velocity_time = omp_get_wtime() - section_start;
  velocity_update_time += step_velocity_time;
  
  // Transfer data back from flat arrays to original data structure
  section_start = omp_get_wtime();
  
  for (int i = 0; i < NumberOfBodies; i++) {
    x[i][0] = pos_x[i];
    x[i][1] = pos_y[i];
    x[i][2] = pos_z[i];
    
    v[i][0] = vel_x[i];
    v[i][1] = vel_y[i];
    v[i][2] = vel_z[i];
  }
  
  step_data_transfer_time += omp_get_wtime() - section_start;
  data_transfer_time += step_data_transfer_time;
  
  // Clean up temporary arrays
  delete[] pos_x;
  delete[] pos_y;
  delete[] pos_z;
  delete[] vel_x;
  delete[] vel_y;
  delete[] vel_z;
  delete[] force_x;
  delete[] force_y;
  delete[] force_z;
  delete[] mass_flat;
  
  t += timeStepSize;
  
  // Calculate total time for this update
  end_time = omp_get_wtime();
  step_total_time = end_time - start_time;
  total_simulation_time += step_total_time;
  
  // Write per-step timing to CSV
  // Only writing on every 10th step to keep file size reasonable
  if (timeStepCounter % 10 == 0) {
    timing_file << std::fixed << std::setprecision(6)
                << timeStepCounter << ","
                << t << ","
                << step_total_time << ","  // This step's time
                << step_force_time << ","
                << step_position_time << ","
                << step_velocity_time << ","
                << step_data_transfer_time << ","
                << (step_force_time/step_total_time*100) << ","
                << (step_position_time/step_total_time*100) << ","
                << (step_velocity_time/step_total_time*100) << ","
                << (step_data_transfer_time/step_total_time*100) << ","
                << NumberOfBodies << std::endl;
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

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
            
  // Also print incremental timing information
  std::cout << "Timing (so far): force calc=" << force_calc_time 
            << "s, position update=" << position_update_time
            << "s, velocity update=" << velocity_update_time 
            << "s, data transfer=" << data_transfer_time << "s"
            << std::endl;
            
  // Write timing data to CSV file
  timing_file << std::fixed << std::setprecision(6)
              << timeStepCounter << ","
              << t << ","
              << total_simulation_time << ","
              << force_calc_time << ","
              << position_update_time << ","
              << velocity_update_time << ","
              << data_transfer_time << ","
              << (force_calc_time/total_simulation_time*100) << ","
              << (position_update_time/total_simulation_time*100) << ","
              << (velocity_update_time/total_simulation_time*100) << ","
              << (data_transfer_time/total_simulation_time*100) << ","
              << NumberOfBodies << std::endl;
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

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}