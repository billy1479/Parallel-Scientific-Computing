#ifndef NBODY_INPUT_READER_H
#define NBODY_INPUT_READER_H

#include <string>
#include <vector>

// A struct to hold simulation parameters
struct SimulationParams {
    double tPlotDelta;
    double tFinal;
    double timeStepSize;
    int numberOfBodies;
    std::vector<double> positions; // x,y,z for each body
    std::vector<double> velocities; // vx,vy,vz for each body
    std::vector<double> masses; // mass for each body
};

class NBodyInputReader {
public:
    // Read simulation parameters from a file
    static SimulationParams readFromFile(const std::string& filename);
    
    // Write simulation parameters to a file
    static void writeToFile(const std::string& filename, const SimulationParams& params);
    
    // Generate parameters using the same logic as create_initial_conditions.py
    static SimulationParams generateInitialConditions(
        int N, int dim, double final_time, double snapshots, double dt,
        double min_mass, double max_mass, const std::string& scenario);
};

#endif // NBODY_INPUT_READER_H