// #include "NBodySimulation.h"
// #include <iostream>
// #include <limits>
// #include <cmath>
// #include <algorithm>
// #include <cstdlib>
// #include <immintrin.h>
// #include <omp.h>
// #include <vector>
// #include <memory>

// class NBodySimulationVectorised : public NBodySimulation {
// public:
//     NBodySimulationVectorised() :
//         t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
//         timeStepCounter(0), timeStepSize(0), maxV(0), minDx(0),
//         snapshotCounter(0) {
//             // Use thread count based on hardware
//             threadCount = omp_get_max_threads();
//             omp_set_num_threads(threadCount);
//         }

//     ~NBodySimulationVectorised() = default; // Use RAII with vectors

//     void setUp(int argc, char** argv) {
//         checkInput(argc, argv);

//         NumberOfBodies = (argc-4) / 7;

//         // Reserve memory to avoid reallocations
//         positions_x.resize(NumberOfBodies);
//         positions_y.resize(NumberOfBodies);
//         positions_z.resize(NumberOfBodies);
//         velocities_x.resize(NumberOfBodies);
//         velocities_y.resize(NumberOfBodies);
//         velocities_z.resize(NumberOfBodies);
//         masses.resize(NumberOfBodies);

//         int readArgument = 1;

//         tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
//         tFinal       = std::stof(argv[readArgument]); readArgument++;
//         timeStepSize = std::stof(argv[readArgument]); readArgument++;

//         for (int i=0; i<NumberOfBodies; i++) {
//             positions_x[i] = std::stof(argv[readArgument]); readArgument++;
//             positions_y[i] = std::stof(argv[readArgument]); readArgument++;
//             positions_z[i] = std::stof(argv[readArgument]); readArgument++;

//             velocities_x[i] = std::stof(argv[readArgument]); readArgument++;
//             velocities_y[i] = std::stof(argv[readArgument]); readArgument++;
//             velocities_z[i] = std::stof(argv[readArgument]); readArgument++;

//             masses[i] = std::stof(argv[readArgument]); readArgument++;

//             if (masses[i]<=0.0) {
//                 std::cerr << "invalid mass for body " << i << std::endl;
//                 exit(-2);
//             }
//         }

//         std::cout << "created setup with " << NumberOfBodies << " bodies"
//                 << " using " << threadCount << " threads" << std::endl;

//         if (tPlotDelta<=0.0) {
//             std::cout << "plotting switched off" << std::endl;
//             tPlot = tFinal + 1.0;
//         }
//         else {
//             std::cout << "plot initial setup plus every " << tPlotDelta
//                     << " time units" << std::endl;
//             tPlot = 0.0;
//         }
//     }

//     void updateBody() {
//         timeStepCounter++;
//         maxV = 0.0;
//         minDx = std::numeric_limits<double>::max();

//         // Use vectors instead of raw arrays
//         std::vector<double> force_x(NumberOfBodies, 0.0);
//         std::vector<double> force_y(NumberOfBodies, 0.0);
//         std::vector<double> force_z(NumberOfBodies, 0.0);
        
//         if (NumberOfBodies == 1) minDx = 0;

//         // Use OpenMP for parallelization of outer loop
//         #pragma omp parallel
//         {
//             // Thread-local variables for reduction
//             double local_minDx = std::numeric_limits<double>::max();

//             // Divide work among threads
//             #pragma omp for schedule(dynamic, 16)
//             for (int i = 0; i < NumberOfBodies; i++) {
//                 // Vectorize inner loop when possible
//                 #pragma omp simd reduction(min:local_minDx)
//                 for (int j = i + 1; j < NumberOfBodies; j++) {
//                     // Calculate forces between bodies i and j
//                     double dx = positions_x[j] - positions_x[i];
//                     double dy = positions_y[j] - positions_y[i];
//                     double dz = positions_z[j] - positions_z[i];
                    
//                     double distSqr = dx*dx + dy*dy + dz*dz;
//                     double dist = std::sqrt(distSqr);
//                     local_minDx = std::min(local_minDx, dist);
                    
//                     double distCubed = dist * distSqr;
//                     double massFactor = masses[i] * masses[j] / distCubed;
                    
//                     double fx = dx * massFactor;
//                     double fy = dy * massFactor;
//                     double fz = dz * massFactor;
                    
//                     // Use atomic updates for thread safety
//                     force_x[i] += fx;
//                     force_y[i] += fy;
//                     force_z[i] += fz;
                    
//                     force_x[j] -= fx;
//                     force_y[j] -= fy;
//                     force_z[j] -= fz;
//                 }
//             }
            
//             // Combine thread-local minDx
//             #pragma omp critical
//             {
//                 minDx = std::min(minDx, local_minDx);
//             }
//         }

//         // Update positions based on current velocities (vectorized)
//         #pragma omp parallel for simd
//         for (int i = 0; i < NumberOfBodies; i++) {
//             positions_x[i] += timeStepSize * velocities_x[i];
//             positions_y[i] += timeStepSize * velocities_y[i];
//             positions_z[i] += timeStepSize * velocities_z[i];
//         }

//         // Update velocities based on forces (vectorized)
//         #pragma omp parallel for simd reduction(max:maxV)
//         for (int i = 0; i < NumberOfBodies; i++) {
//             velocities_x[i] += timeStepSize * force_x[i] / masses[i];
//             velocities_y[i] += timeStepSize * force_y[i] / masses[i];
//             velocities_z[i] += timeStepSize * force_z[i] / masses[i];
            
//             // Calculate velocity magnitude for statistics
//             double v_squared = velocities_x[i]*velocities_x[i] + 
//                               velocities_y[i]*velocities_y[i] + 
//                               velocities_z[i]*velocities_z[i];
//             double v_mag = std::sqrt(v_squared);
//             maxV = std::max(maxV, v_mag);
//         }

//         t += timeStepSize;
//     }

//     bool hasReachedEnd() {
//         return t > tFinal;
//     }

//     void takeSnapshot() {
//         if (t >= tPlot) {
//             printParaviewSnapshot();
//             printSnapshotSummary();
//             tPlot += tPlotDelta;
//         }
//     }
    
//     void openParaviewVideoFile() {
//         videoFile.open("paraview-output/result.pvd");
//         videoFile << "<?xml version=\"1.0\"?>" << std::endl
//                 << "<VTKFile type=\"Collection\""
//                 << " version=\"0.1\""
//                 << " byte_order=\"LittleEndian\""
//                 << " compressor=\"vtkZLibDataCompressor\">" << std::endl
//                 << "<Collection>";
//     }
    
//     void closeParaviewVideoFile() {
//         videoFile << "</Collection>"
//                 << "</VTKFile>" << std::endl;
//         videoFile.close();
//     }

//     void printParaviewSnapshot() {
//         static int counter = -1;
//         counter++;
//         std::stringstream filename, filename_nofolder;
//         filename << "paraview-output/result-" << counter << ".vtp";
//         filename_nofolder << "result-" << counter << ".vtp";
//         std::ofstream out(filename.str().c_str());
        
//         out << "<VTKFile type=\"PolyData\" >" << std::endl
//             << "<PolyData>" << std::endl
//             << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
//             << "  <Points>" << std::endl
//             << "   <DataArray type=\"Float64\""
//             << " NumberOfComponents=\"3\""
//             << " format=\"ascii\">";

//         // Write positions using separate coordinate arrays
//         for (int i = 0; i < NumberOfBodies; i++) {
//             out << positions_x[i] << " "
//                 << positions_y[i] << " "
//                 << positions_z[i] << " ";
//         }

//         out << "   </DataArray>" << std::endl
//             << "  </Points>" << std::endl
//             << " </Piece>" << std::endl
//             << "</PolyData>" << std::endl
//             << "</VTKFile>" << std::endl;

//         out.close();

//         videoFile << "<DataSet timestep=\"" << counter
//                 << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
//                 << "\"/>" << std::endl;
//     }

//     void printSnapshotSummary() {
//         std::cout << "plot next snapshot"
//                 << ",\t time step=" << timeStepCounter
//                 << ",\t t=" << t
//                 << ",\t dt=" << timeStepSize
//                 << ",\t v_max=" << maxV
//                 << ",\t dx_min=" << minDx
//                 << std::endl;
//     }

//     void printSummary() {
//         std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
//         std::cout << "Position of first remaining object: "
//                 << positions_x[0] << ", " << positions_y[0] << ", " << positions_z[0] << std::endl;
//     }

// private:
//     // Use STL vectors instead of raw pointers for automatic memory management
//     std::vector<double> positions_x;
//     std::vector<double> positions_y;
//     std::vector<double> positions_z;
//     std::vector<double> velocities_x;
//     std::vector<double> velocities_y;
//     std::vector<double> velocities_z;
//     std::vector<double> masses;

//     double t;
//     double tFinal;
//     double tPlot;
//     double tPlotDelta;
//     int NumberOfBodies;
//     int timeStepCounter;
//     double timeStepSize;
//     int threadCount;

//     // Statistics and monitoring
//     double maxV;
//     double minDx;
    
//     // Output handling
//     std::ofstream videoFile;
//     int snapshotCounter;
// };


#include "NBodySimulation.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>

// Memory-aligned structure for better cache performance
struct alignas(64) Body {
    double x, y, z;        // Position
    double vx, vy, vz;     // Velocity
    double mass;           // Mass
    double padding;        // Maintain 64-byte alignment
};

class NBodySimulationVectorised : public NBodySimulation {
public:
    NBodySimulationVectorised() :
        t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
        timeStepCounter(0), timeStepSize(0), maxV(0), minDx(0),
        snapshotCounter(0) {
            bodies = nullptr;
        }

    ~NBodySimulationVectorised() {
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