#include "NBodySimulation.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <immintrin.h>
#include <omp.h>

class NBodySimulationVectorised : public NBodySimulation {
    public:
        NBodySimulationVectorised() :
            t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
            timeStepCounter(0), timeStepSize(0),
            x_pos(nullptr), y_pos(nullptr), z_pos(nullptr),
            x_vel(nullptr), y_vel(nullptr), z_vel(nullptr),
            mass(nullptr), maxV(0), minDx(0),
            snapshotCounter(0) {}

        ~NBodySimulationVectorised() {
            delete[] x_pos;
            delete[] y_pos;
            delete[] z_pos;
            delete[] x_vel;
            delete[] y_vel;
            delete[] z_vel;
            delete[] mass;
        }

        void setUp(int argc, char** argv) {
            checkInput(argc, argv);

            NumberOfBodies = (argc-4) / 7;

            // Allocate 1D arrays
            x_pos = new double[NumberOfBodies];
            y_pos = new double[NumberOfBodies];
            z_pos = new double[NumberOfBodies];
            x_vel = new double[NumberOfBodies];
            y_vel = new double[NumberOfBodies];
            z_vel = new double[NumberOfBodies];
            mass = new double[NumberOfBodies];

            int readArgument = 1;

            tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
            tFinal       = std::stof(argv[readArgument]); readArgument++;
            timeStepSize = std::stof(argv[readArgument]); readArgument++;

            for (int i=0; i<NumberOfBodies; i++) {
                x_pos[i] = std::stof(argv[readArgument]); readArgument++;
                y_pos[i] = std::stof(argv[readArgument]); readArgument++;
                z_pos[i] = std::stof(argv[readArgument]); readArgument++;

                x_vel[i] = std::stof(argv[readArgument]); readArgument++;
                y_vel[i] = std::stof(argv[readArgument]); readArgument++;
                z_vel[i] = std::stof(argv[readArgument]); readArgument++;

                mass[i] = std::stof(argv[readArgument]); readArgument++;

                if (mass[i]<=0.0) {
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


        #pragma omp declare simd // SIMD enabled version of the square root function
        double simd_sqrt(float x) {
            return sqrt(x);
        }

        inline double simd_max(double a, double b) { // Open mp max function
            __m128d va = _mm_set_sd(a);
            __m128d vb = _mm_set_sd(b);
            __m128d vmax = _mm_max_sd(va, vb);
            return _mm_cvtsd_f64(vmax);
        }

        inline double simd_min(double a, double b) { // Open mp min function
            __m128d va = _mm_set_sd(a);
            __m128d vb = _mm_set_sd(b);
            __m128d vmin = _mm_min_sd(va, vb);
            return _mm_cvtsd_f64(vmin);
        }

        void updateBody () {
            timeStepCounter++;
            maxV   = 0.0;
            minDx  = std::numeric_limits<double>::max();

            double* force0 = new double[NumberOfBodies]();
            double* force1 = new double[NumberOfBodies]();
            double* force2 = new double[NumberOfBodies]();

            if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

            #pragma omp simd reduction(+: force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies]) // HAVE OPTIMISED BY USING SIMD TO vectorise THE LOOP
            for (int i=0; i<NumberOfBodies; i++) { // HAVE OPTIMISED BY USING VARIABLES TO STORE THE FORCE CALCULATION
                for (int j=i+1; j<NumberOfBodies; j++) {
                    if(i!=j){
                        double f_0 = force_calculation(i,j,0);
                        double f_1 = force_calculation(i,j,1);
                        double f_2 = force_calculation(i,j,2);
                        // x,y,z forces acting on particle i.
                        force0[i] += f_0;
                        force1[i] += f_1;
                        force2[i] += f_2;
                        // x,y,z symmetric forces acting on particle j.
                        force0[j] -= f_0;
                        force1[j] -= f_1;
                        force2[j] -= f_2;

                    }
                }
            }

            for (int i=0; i < NumberOfBodies; i++){
                x_pos[i] += timeStepSize * x_vel[i];
                y_pos[i] += timeStepSize * y_vel[i];
                z_pos[i] += timeStepSize * z_vel[i];
            }


           for (int i=0; i < NumberOfBodies; i++) {
                x_vel[i] += timeStepSize * force0[i] / mass[i];
                y_vel[i] += timeStepSize * force1[i] / mass[i];
                z_vel[i] += timeStepSize * force2[i] / mass[i];

                maxV = std::max(maxV, std::sqrt(x_vel[i]*x_vel[i] + y_vel[i]*y_vel[i] + z_vel[i]*z_vel[i]));
            }

            t += timeStepSize;

            delete[] force0;
            delete[] force1;
            delete[] force2;
        }

        #pragma omp declare simd
        double force_calculation (int i, int j, int direction){
            const double dx = x_pos[i] - x_pos[j];
            const double dy = y_pos[i] - y_pos[j];
            const double dz = z_pos[i] - z_pos[j];
            
            const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);

            minDx = std::min(minDx,distance);

            double coord_diff;
            switch(direction) {
                case 0: coord_diff = dx; break;
                case 1: coord_diff = dy; break;
                case 2: coord_diff = dz; break;
                default: coord_diff = 0;
            }

            return coord_diff * mass[i] * mass[j] / (distance * distance * distance);
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

            // Write positions using separate coordinate arrays
            for (int i = 0; i < NumberOfBodies; i++) {
                out << x_pos[i] << " "
                    << y_pos[i] << " "
                    << z_pos[i] << " ";
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
                    << x_pos[0] << ", " << y_pos[0] << ", " << z_pos[0] << std::endl;
        }

    private:
        double* x_pos;
        double* y_pos;
        double* z_pos;
        double* x_vel;
        double* y_vel;
        double* z_vel;

        double* mass;

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
