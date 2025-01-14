#include <iomanip>
#include <immintrin.h>
#include <omp.h>


#include "NBodySimulationVectorised.cpp"

/**
 * You can compile this file with
 *   make step-2-g++   // Uses the GNU Compiler Collection.
 *   make step-2-icpx  // Uses the Intel compiler.
 * and run it with
 *   ./step-2-g++
 *   ./step-2-icpx
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

class NBodySimulationParallelised : public NBodySimulationVectorised {
    public:
        NBodySimulationParallelised() :
            t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
            timeStepCounter(0), timeStepSize(0),
            x_pos(nullptr), y_pos(nullptr), z_pos(nullptr),
            x_vel(nullptr), y_vel(nullptr), z_vel(nullptr),
            mass(nullptr), maxV(0), minDx(0),
            snapshotCounter(0) {
            // Set number of threads at construction
            num_threads = omp_get_max_threads();
            omp_set_num_threads(num_threads);
        }

        ~NBodySimulationParallelised() {
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
            x_pos = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            y_pos = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            z_pos = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            x_vel = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            y_vel = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            z_vel = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
            mass = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));

            int readArgument = 1;

            tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
            tFinal       = std::stof(argv[readArgument]); readArgument++;
            timeStepSize = std::stof(argv[readArgument]); readArgument++;

            for (int i=0; i<NumberOfBodies; i++) {
                x_pos[i] = std::stof(argv[readArgument + i*7]);
                y_pos[i] = std::stof(argv[readArgument + i*7 + 1]);
                z_pos[i] = std::stof(argv[readArgument + i*7 + 2]);
                x_vel[i] = std::stof(argv[readArgument + i*7 + 3]);
                y_vel[i] = std::stof(argv[readArgument + i*7 + 4]);
                z_vel[i] = std::stof(argv[readArgument + i*7 + 5]);
                mass[i] = std::stof(argv[readArgument + i*7 + 6]);

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
        void updateBody () {

          timeStepCounter++;
          maxV   = 0.0;
          minDx  = std::numeric_limits<double>::max();

          // force0 = force along x direction
          // force1 = force along y direction
          // force2 = force along z direction
          double* force0 = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
          double* force1 = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));
          double* force2 = static_cast<double*>(_mm_malloc(NumberOfBodies * sizeof(double), 32));

          if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

          // Set number of threads
          omp_set_num_threads(5);

          // #pragma omp parallel for schedule(dynamic) reduction(+: force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies])
          for (int i = 0; i < NumberOfBodies; i++) {
            for (int j = i + 1; j < NumberOfBodies; j++) {
              if (i != j) {
                double fx = force_calculation(i, j, 0);
                double fy = force_calculation(i, j, 1);
                double fz = force_calculation(i, j, 2);

                // #pragma omp atomic
                force0[i] += fx;
                // #pragma omp atomic
                force1[i] += fy;
                // #pragma omp atomic
                force2[i] += fz;

                // #pragma omp atomic
                force0[j] -= fx;
                // #pragma omp atomic
                force1[j] -= fy;
                // #pragma omp atomic
                force2[j] -= fz;
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

          _mm_free(force0);
          _mm_free(force1);
          _mm_free(force2);
        }
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
      int num_threads;

      double maxV;
      double minDx;
      
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
  NBodySimulationParallelised nbs;
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
