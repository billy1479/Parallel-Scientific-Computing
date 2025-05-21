#include <iostream>
#include <iomanip>
#include <string>
#include "NBodySimulation.h"
#include "NBodyInputReader.h"

void printUsage(const char* programName) {
  std::cout << "Usage: " << programName << " [options]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  --input FILENAME   Read simulation parameters from input file" << std::endl;
  std::cout << "  or use original command line parameters:" << std::endl;
  std::cout << "    plot-time final-time dt objects" << std::endl;
}

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      printUsage(argv[0]);
      return 1;
    }

    std::cout << std::setprecision(15);

    NBodySimulation nbs;

    // Check for --input argument
    if (std::string(argv[1]) == "--input") {
      if (argc < 3) {
        std::cerr << "Error: Missing input file after --input" << std::endl;
        return 1;
      }

      try {
        // Read parameters from file
        SimulationParams params = NBodyInputReader::readFromFile(argv[2]);
        
        // Set up simulation using params
        nbs.setUpFromParams(params);
        
        std::cout << "Successfully loaded " << params.numberOfBodies 
                  << " bodies from input file: " << argv[2] << std::endl;
      }
      catch (const std::exception& e) {
        std::cerr << "Error reading input file: " << e.what() << std::endl;
        return 1;
      }
    } else {
      // Use command-line arguments
      nbs.setUp(argc, argv);
    }

    nbs.openParaviewVideoFile();
    nbs.takeSnapshot();

    while (!nbs.hasReachedEnd()) {
      nbs.updateBody();
      nbs.takeSnapshot();
    }

    nbs.printSummary();
    nbs.closeParaviewVideoFile();
  } 
  catch (const std::runtime_error& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
  catch (int e) {
    // Handle legacy error code
    return e;
  }

  return 0;
}