#include "NBodyInputReader.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <random>
#include <stdexcept>

SimulationParams NBodyInputReader::readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    
    SimulationParams params;
    std::string line;
    
    // Read header line with basic parameters
    std::getline(file, line);
    std::istringstream header(line);
    header >> params.tPlotDelta >> params.tFinal >> params.timeStepSize >> params.numberOfBodies;
    
    // Resize vectors for positions, velocities, and masses
    params.positions.resize(params.numberOfBodies * 3);  // x,y,z for each body
    params.velocities.resize(params.numberOfBodies * 3); // vx,vy,vz for each body
    params.masses.resize(params.numberOfBodies);        // mass for each body
    
    // Read body data
    for (int i = 0; i < params.numberOfBodies; i++) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("Error reading body data at index " + std::to_string(i));
        }
        
        std::istringstream bodyData(line);
        
        // Read position x,y,z
        bodyData >> params.positions[i*3] >> params.positions[i*3+1] >> params.positions[i*3+2];
        
        // Read velocity vx,vy,vz
        bodyData >> params.velocities[i*3] >> params.velocities[i*3+1] >> params.velocities[i*3+2];
        
        // Read mass
        bodyData >> params.masses[i];
        
        if (params.masses[i] <= 0.0) {
            throw std::runtime_error("Invalid mass for body " + std::to_string(i));
        }
    }
    
    return params;
}

void NBodyInputReader::writeToFile(const std::string& filename, const SimulationParams& params) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    
    // Write header with basic parameters
    file << params.tPlotDelta << " " << params.tFinal << " "
         << params.timeStepSize << " " << params.numberOfBodies << std::endl;
    
    // Write body data
    for (int i = 0; i < params.numberOfBodies; i++) {
        // Write position x,y,z
        file << params.positions[i*3] << " " 
             << params.positions[i*3+1] << " "
             << params.positions[i*3+2] << " ";
        
        // Write velocity vx,vy,vz
        file << params.velocities[i*3] << " "
             << params.velocities[i*3+1] << " "
             << params.velocities[i*3+2] << " ";
        
        // Write mass
        file << params.masses[i] << std::endl;
    }
}

SimulationParams NBodyInputReader::generateInitialConditions(
    int N, int dim, double final_time, double snapshots, double dt,
    double min_mass, double max_mass, const std::string& scenario) {
    
    SimulationParams params;
    params.tPlotDelta = snapshots;
    params.tFinal = final_time;
    params.timeStepSize = dt;
    
    // Determine dimensions
    int N_x = N;
    int N_y = (dim >= 2) ? N : 1;
    int N_z = (dim >= 3) ? N : 1;
    
    // Calculate number of bodies
    params.numberOfBodies = N_x * N_y * N_z;
    
    // Resize arrays
    params.positions.resize(params.numberOfBodies * 3);
    params.velocities.resize(params.numberOfBodies * 3);
    params.masses.resize(params.numberOfBodies);
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_noise(-0.45, 0.45);
    std::uniform_real_distribution<> mass_dist(min_mass, max_mass);
    
    // Generate body data
    int idx = 0;
    for (int x = 0; x < N_x; x++) {
        for (int y = 0; y < N_y; y++) {
            for (int z = 0; z < N_z; z++) {
                double xPos, yPos, zPos;
                double xVel = 0.0, yVel = 0.0, zVel = 0.0;
                
                if (scenario == "random-grid") {
                    xPos = (pos_noise(gen)) * 0.9 * (1.0/N_x) + x * (1.0/N_x);
                    yPos = (pos_noise(gen)) * 0.9 * (1.0/N_y) + y * (1.0/N_y);
                    zPos = (pos_noise(gen)) * 0.9 * (1.0/N_z) + z * (1.0/N_z);
                } else if (scenario == "no-noise") {
                    xPos = x * (1.0/N_x);
                    yPos = y * (1.0/N_y);
                    zPos = z * (1.0/N_z);
                } else if (scenario == "shock") {
                    xPos = (pos_noise(gen)) * 0.9 * (1.0/N_x) + x * (1.0/N_x);
                    yPos = (pos_noise(gen)) * 0.9 * (1.0/N_y) + y * (1.0/N_y);
                    zPos = (pos_noise(gen)) * 0.9 * (1.0/N_z) + z * (1.0/N_z);
                    
                    double dist = std::sqrt((xPos-0.5)*(xPos-0.5) + 
                                           (yPos-0.5)*(yPos-0.5) + 
                                           (zPos-0.5)*(zPos-0.5));
                    if (dist < 0.1) {
                        xVel = (xPos-0.5) / (dist+0.00001);
                        yVel = (yPos-0.5) / (dist+0.00001);
                        zVel = (zPos-0.5) / (dist+0.00001);
                    }
                } else {
                    throw std::runtime_error("Unknown scenario: " + scenario);
                }
                
                double mass = mass_dist(gen);
                
                // Store data
                params.positions[idx*3] = xPos;
                params.positions[idx*3+1] = yPos;
                params.positions[idx*3+2] = zPos;
                
                params.velocities[idx*3] = xVel;
                params.velocities[idx*3+1] = yVel;
                params.velocities[idx*3+2] = zVel;
                
                params.masses[idx] = mass;
                
                idx++;
            }
        }
    }
    
    return params;
}