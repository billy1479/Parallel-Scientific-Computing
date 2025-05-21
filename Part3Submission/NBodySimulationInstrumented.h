#ifndef NBODY_SIMULATION_INSTRUMENTED_H
#define NBODY_SIMULATION_INSTRUMENTED_H

#include "NBodySimulation.h"
#include <omp.h>
#include <fstream>

/**
 * Extended N-body simulation class with performance instrumentation.
 * This class inherits from NBodySimulation and adds timing measurements
 * for performance analysis.
 */
class NBodySimulationInstrumented : public NBodySimulation {
private:
    // Timing variables for the performance analysis
    double force_calc_time;
    double position_update_time;
    double velocity_update_time;
    double total_simulation_time;

    // Per-step timing variables
    double step_total_time;
    double step_force_time;
    double step_position_time;
    double step_velocity_time;

    // File for logging timing data
    std::ofstream timing_file;

public:
    /**
     * Constructor that initializes timing variables and creates the timing CSV file.
     */
    NBodySimulationInstrumented();
    
    /**
     * Destructor that prints timing summary and closes the timing file.
     */
    ~NBodySimulationInstrumented();
    
    /**
     * Override methods that need instrumentation
     */
    void setUp(int argc, char** argv);
    double force_calculation(int i, int j, int direction);
    void updateBody();
    void printSnapshotSummary();
    
    /**
     * Passthrough methods that call the base class
     */
    void openParaviewVideoFile();
    void takeSnapshot();
    bool hasReachedEnd();
    void printSummary();
    void closeParaviewVideoFile();
};

#endif // NBODY_SIMULATION_INSTRUMENTED_H