#include "NBodySimulation.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>
#include <iomanip>

/**
 * You can compile this file with
 *   make step-3.
 * and run it with
 *   ./step-3
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

 // Memory-aligned structure for better cache performance
struct alignas(64) Body {
  double x, y, z;        // Position
  double vx, vy, vz;     // Velocity
  double mass;           // Mass
  double padding;        // Maintain 64-byte alignment
};

class NBodySimulationMolecularForces : public NBodySimulation {
  public:
  NBodySimulationMolecularForces() :
      t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
      timeStepCounter(0), timeStepSize(0), maxV(0), minDx(0),
      snapshotCounter(0) {
          bodies = nullptr;
      }

  ~NBodySimulationMolecularForces() {
      if (bodies) {
          _mm_free(bodies);  // Free aligned memory
      }
  }

  void setUp(int argc, char** argv) {
      checkInput(argc, argv);

      NumberOfBodies = (argc-4) / 7;

      // Allocate aligned memory for bodies
      bodies = static_cast<Body*>(_mm_malloc(NumberOfBodies * sizeof(Body), 64));
      if (!bodies) {
          std::cerr << "Failed to allocate aligned memory" << std::endl;
          exit(-1);
      }

      int readArgument = 1;

      tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
      tFinal       = std::stof(argv[readArgument]); readArgument++;
      timeStepSize = std::stof(argv[readArgument]); readArgument++;

      for (int i=0; i<NumberOfBodies; i++) {
          bodies[i].x = std::stof(argv[readArgument]); readArgument++;
          bodies[i].y = std::stof(argv[readArgument]); readArgument++;
          bodies[i].z = std::stof(argv[readArgument]); readArgument++;

          bodies[i].vx = std::stof(argv[readArgument]); readArgument++;
          bodies[i].vy = std::stof(argv[readArgument]); readArgument++;
          bodies[i].vz = std::stof(argv[readArgument]); readArgument++;

          bodies[i].mass = std::stof(argv[readArgument]); readArgument++;

          if (bodies[i].mass <= 0.0) {
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

  void updateBody() {
      timeStepCounter++;
      maxV   = 0.0;
      minDx  = std::numeric_limits<double>::max();

      // Allocate force arrays
      double* force0 = new double[NumberOfBodies]();
      double* force1 = new double[NumberOfBodies]();
      double* force2 = new double[NumberOfBodies]();

      if (NumberOfBodies == 1) minDx = 0;

      for (int i = 0; i < NumberOfBodies; i++) {
          for (int j = i + 1; j < NumberOfBodies; j++) {
              if (i != j) {
                  double f_0 = force_calculation(i, j, 0);
                  double f_1 = force_calculation(i, j, 1);
                  double f_2 = force_calculation(i, j, 2);

                  force0[i] += f_0;
                  force1[i] += f_1;
                  force2[i] += f_2;

                  force0[j] -= f_0;
                  force1[j] -= f_1;
                  force2[j] -= f_2;
              }
          }
      }

      for (int i = 0; i < NumberOfBodies; i++) {
          bodies[i].x += timeStepSize * bodies[i].vx;
          bodies[i].y += timeStepSize * bodies[i].vy;
          bodies[i].z += timeStepSize * bodies[i].vz;

          bodies[i].vx += timeStepSize * force0[i] / bodies[i].mass;
          bodies[i].vy += timeStepSize * force1[i] / bodies[i].mass;
          bodies[i].vz += timeStepSize * force2[i] / bodies[i].mass;

          maxV = std::max(maxV, std::sqrt(bodies[i].vx*bodies[i].vx + 
                                        bodies[i].vy*bodies[i].vy + 
                                        bodies[i].vz*bodies[i].vz));
      }

      t += timeStepSize;

      delete[] force0;
      delete[] force1;
      delete[] force2;
  }

  double force_calculation(int i, int j, int direction) {
      double dx = bodies[j].x - bodies[i].x;
      double dy = bodies[j].y - bodies[i].y;
      double dz = bodies[j].z - bodies[i].z;
      
      const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
      const double distance3 = distance * distance * distance;
      minDx = std::min(minDx, distance);

      double coord_diff;
      switch(direction) {
          case 0: coord_diff = dx; break;
          case 1: coord_diff = dy; break;
          case 2: coord_diff = dz; break;
          default: coord_diff = 0;
      }

      return coord_diff * bodies[i].mass * bodies[j].mass / distance3;
  }

  bool hasReachedEnd() {
      return t > tFinal;
  }

  void takeSnapshot() {
      if (t >= tPlot) {
          printParaviewSnapshot();
          printSnapshotSummary();
          tPlot += tPlotDelta;
      }
  }
  
  void openParaviewVideoFile() {
      videoFile.open("paraview-output/result.pvd");
      videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\""
              << " version=\"0.1\""
              << " byte_order=\"LittleEndian\""
              << " compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
  }
  
  void closeParaviewVideoFile() {
      videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
      videoFile.close();
  }

  void printParaviewSnapshot() {
      static int counter = -1;
      counter++;
      std::stringstream filename, filename_nofolder;
      filename << "paraview-output/result-" << counter << ".vtp";
      filename_nofolder << "result-" << counter << ".vtp";
      std::ofstream out(filename.str().c_str());
      
      out << "<VTKFile type=\"PolyData\" >" << std::endl
          << "<PolyData>" << std::endl
          << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
          << "  <Points>" << std::endl
          << "   <DataArray type=\"Float64\""
          << " NumberOfComponents=\"3\""
          << " format=\"ascii\">";

      for (int i = 0; i < NumberOfBodies; i++) {
          out << bodies[i].x << " "
              << bodies[i].y << " "
              << bodies[i].z << " ";
      }

      out << "   </DataArray>" << std::endl
          << "  </Points>" << std::endl
          << " </Piece>" << std::endl
          << "</PolyData>" << std::endl
          << "</VTKFile>" << std::endl;

      out.close();

      videoFile << "<DataSet timestep=\"" << counter
              << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
              << "\"/>" << std::endl;
  }

  void printSnapshotSummary() {
      std::cout << "plot next snapshot"
              << ",\t time step=" << timeStepCounter
              << ",\t t=" << t
              << ",\t dt=" << timeStepSize
              << ",\t v_max=" << maxV
              << ",\t dx_min=" << minDx
              << std::endl;
  }

  void printSummary() {
      std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
      std::cout << "Position of first remaining object: "
              << bodies[0].x << ", " << bodies[0].y << ", " << bodies[0].z << std::endl;
  }

private:
  Body* bodies;           // Memory-aligned array of bodies

  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;
  int NumberOfBodies;
  int timeStepCounter;
  double timeStepSize;

  // Statistics and monitoring
  double maxV;
  double minDx;
  
  // Output handling
  std::ofstream videoFile;
  int snapshotCounter;
};

/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */
int main (int argc, char** argv) {

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationMolecularForces nbs;
  nbs.setUp(argc,argv);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
