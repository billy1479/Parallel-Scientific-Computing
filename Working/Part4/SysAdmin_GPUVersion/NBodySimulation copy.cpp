#include "NBodySimulation.h"
#include <omp.h>  // Add OpenMP header
#include <iomanip> // For std::setprecision

// Add timing variables to the class
double force_calc_time = 0.0;
double position_update_time = 0.0;
double velocity_update_time = 0.0;
double total_simulation_time = 0.0;

NBodySimulation::NBodySimulation () :
  force_calc_time(0.0), position_update_time(0.0), velocity_update_time(0.0), total_simulation_time(0.0),
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  x(nullptr), v(nullptr), mass(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {
  
  // Create and initialize the timing CSV file
  timing_file.open("nbody_timing.csv");
  // Write CSV header
  timing_file << "TimeStep,SimulationTime,TotalElapsedTime,ForceCalcTime,PositionUpdateTime,VelocityUpdateTime,ForceCalcPercentage,PositionUpdatePercentage,VelocityUpdatePercentage,NumberOfBodies" << std::endl;
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
  
  // Write final summary to CSV and close file
  timing_file << "SUMMARY,,"
              << total_simulation_time << ","
              << force_calc_time << ","
              << position_update_time << ","
              << velocity_update_time << ","
              << (force_calc_time/total_simulation_time*100) << ","
              << (position_update_time/total_simulation_time*100) << ","
              << (velocity_update_time/total_simulation_time*100) << ","
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
  double start_time, end_time, section_start;
  double step_force_time = 0.0;
  double step_position_time = 0.0;
  double step_velocity_time = 0.0;
  double step_total_time = 0.0;
  
  // Start measuring total time for this update
  start_time = omp_get_wtime();

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

  // Force calculation section - typically the most expensive part
  section_start = omp_get_wtime();
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
  step_force_time = omp_get_wtime() - section_start;
  force_calc_time += step_force_time;
  
  // Position update section
  section_start = omp_get_wtime();
  for (int i=0; i < NumberOfBodies; i++){
    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];
  }
  step_position_time = omp_get_wtime() - section_start;
  position_update_time += step_position_time;
  
  // Velocity update section
  section_start = omp_get_wtime();
  for (int i=0; i < NumberOfBodies; i++){
    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

    maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
  }
  step_velocity_time = omp_get_wtime() - section_start;
  velocity_update_time += step_velocity_time;
  
  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
  
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
                << (step_force_time/step_total_time*100) << ","
                << (step_position_time/step_total_time*100) << ","
                << (step_velocity_time/step_total_time*100) << ","
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
            
  // Also print incremental timing information
  std::cout << "Timing (so far): force calc=" << force_calc_time 
            << "s, position update=" << position_update_time
            << "s, velocity update=" << velocity_update_time << "s"
            << std::endl;
            
  // Write timing data to CSV file
  // Format: TimeStep, SimulationTime, TotalElapsedTime, ForceCalcTime, PositionUpdateTime, VelocityUpdateTime, 
  //         ForceCalcPercentage, PositionUpdatePercentage, VelocityUpdatePercentage, NumberOfBodies
  timing_file << std::fixed << std::setprecision(6)
              << timeStepCounter << ","
              << t << ","
              << total_simulation_time << ","
              << force_calc_time << ","
              << position_update_time << ","
              << velocity_update_time << ","
              << (force_calc_time/total_simulation_time*100) << ","
              << (position_update_time/total_simulation_time*100) << ","
              << (velocity_update_time/total_simulation_time*100) << ","
              << NumberOfBodies << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;
}