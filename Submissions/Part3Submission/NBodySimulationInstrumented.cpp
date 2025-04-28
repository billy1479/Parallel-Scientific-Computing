#include "NBodySimulationInstrumented.h"
#include <iomanip> // For std::setprecision

NBodySimulationInstrumented::NBodySimulationInstrumented() :
    force_calc_time(0.0), position_update_time(0.0), velocity_update_time(0.0), total_simulation_time(0.0),
    step_total_time(0.0), step_force_time(0.0), step_position_time(0.0), step_velocity_time(0.0) {
    // Create and initialize the timing CSV file
    timing_file.open("nbody_timing.csv");
    timing_file << "TimeStep,SimulationTime,TotalElapsedTime,ForceCalcTime,PositionUpdateTime,VelocityUpdateTime,ForceCalcPercentage,PositionUpdatePercentage,VelocityUpdatePercentage,NumberOfBodies" << std::endl;
}

void NBodySimulationInstrumented::setUp(int argc, char** argv) {
  NBodySimulation::setUp(argc, argv);
}

void NBodySimulationInstrumented::openParaviewVideoFile() {
  NBodySimulation::openParaviewVideoFile();
}

void NBodySimulationInstrumented::takeSnapshot() {
  NBodySimulation::takeSnapshot();
}

bool NBodySimulationInstrumented::hasReachedEnd() {
  return NBodySimulation::hasReachedEnd();
}

void NBodySimulationInstrumented::printSummary() {
  NBodySimulation::printSummary();
}

void NBodySimulationInstrumented::closeParaviewVideoFile() {
  NBodySimulation::closeParaviewVideoFile();
}

double NBodySimulationInstrumented::force_calculation(int i, int j, int direction) {
    // Call the base class implementation first
    return NBodySimulation::force_calculation(i, j, direction);
}

NBodySimulationInstrumented::~NBodySimulationInstrumented() {
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

void NBodySimulationInstrumented::updateBody() {
    double start_time, end_time, section_start;
    
    // Start measuring total time for this update
    start_time = omp_get_wtime();
  
    timeStepCounter++;
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();
  
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
  
      maxV = std::max(maxV, std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]));
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

void NBodySimulationInstrumented::printSnapshotSummary() {
    // Call the base class implementation first
    NBodySimulation::printSnapshotSummary();

    // Add the instrumentation-specific output
    std::cout << "Timing (so far): force calc=" << force_calc_time 
              << "s, position update=" << position_update_time
              << "s, velocity update=" << velocity_update_time << "s"
              << std::endl;

    // Write timing data to CSV file
    timing_file << std::fixed << std::setprecision(6)
                << timeStepCounter << ","
                << t << ","
                << step_total_time << ","
                << step_force_time << ","
                << step_position_time << ","
                << step_velocity_time << ","
                << (step_force_time / step_total_time * 100) << ","
                << (step_position_time / step_total_time * 100) << ","
                << (step_velocity_time / step_total_time * 100) << ","
                << NumberOfBodies << std::endl;
}