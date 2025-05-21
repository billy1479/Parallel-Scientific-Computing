#include "NBodyInputReader.h"
#include <iostream>
#include <string>

void printUsage() {
    std::cout << "Usage: ./initial_conditions_generator [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  --final-time X       Simulation runs from 0 through final time X" << std::endl;
    std::cout << "  --snapshots X        Simulation writes a snapshot every X time units" << std::endl;
    std::cout << "  --dt X               Time step size" << std::endl;
    std::cout << "  --min-mass X         Minimal mass of particles" << std::endl;
    std::cout << "  --max-mass X         Maximal mass of particles" << std::endl;
    std::cout << "  --dim X              Dimensions (1, 2, or 3)" << std::endl;
    std::cout << "  --N X                Grid size (N x N x N particles)" << std::endl;
    std::cout << "  --scenario X         Scenario (random-grid, shock, or no-noise)" << std::endl;
    std::cout << "  --output X           Output file name" << std::endl;
}

int main(int argc, char** argv) {
    if (argc <= 1) {
        printUsage();
        return 1;
    }
    
    // Default values
    double final_time = 10.0;
    double snapshots = 0.1;
    double dt = 0.001;
    double min_mass = 0.1;
    double max_mass = 1.0;
    int dim = 3;
    int N = 10;
    std::string scenario = "random-grid";
    std::string output_file = "initial_conditions.nbody";
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--final-time" && i + 1 < argc) {
            final_time = std::stod(argv[++i]);
        } else if (arg == "--snapshots" && i + 1 < argc) {
            snapshots = std::stod(argv[++i]);
        } else if (arg == "--dt" && i + 1 < argc) {
            dt = std::stod(argv[++i]);
        } else if (arg == "--min-mass" && i + 1 < argc) {
            min_mass = std::stod(argv[++i]);
        } else if (arg == "--max-mass" && i + 1 < argc) {
            max_mass = std::stod(argv[++i]);
        } else if (arg == "--dim" && i + 1 < argc) {
            dim = std::stoi(argv[++i]);
            if (dim < 1 || dim > 3) {
                std::cerr << "Error: Dimension must be 1, 2, or 3" << std::endl;
                return 1;
            }
        } else if (arg == "--N" && i + 1 < argc) {
            N = std::stoi(argv[++i]);
        } else if (arg == "--scenario" && i + 1 < argc) {
            scenario = argv[++i];
            if (scenario != "random-grid" && scenario != "shock" && scenario != "no-noise") {
                std::cerr << "Error: Unknown scenario. Supported values are random-grid, shock, and no-noise" << std::endl;
                return 1;
            }
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            printUsage();
            return 1;
        }
    }
    
    try {
        // Generate initial conditions
        SimulationParams params = NBodyInputReader::generateInitialConditions(
            N, dim, final_time, snapshots, dt, min_mass, max_mass, scenario);
        
        // Write to file
        NBodyInputReader::writeToFile(output_file, params);
        
        int totalBodies = params.numberOfBodies;
        std::cout << "Created setup with " << totalBodies << " bodies" << std::endl;
        std::cout << "Wrote to file: " << output_file << std::endl;
        std::cout << "To run the simulation, use: ./step-5-gpu --input " << output_file << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}