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

          // force0 = force along x direction
          // force1 = force along y direction
          // force2 = force along z direction
          double* force0 = new double[NumberOfBodies]();
          double* force1 = new double[NumberOfBodies]();
          double* force2 = new double[NumberOfBodies]();

          if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

          // Set number of threads
          omp_set_num_threads(5);

          #pragma omp parallel for schedule(dynamic) reduction(+: force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies])
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

          #pragma omp parallel for
          for (int i = 0; i < NumberOfBodies; i++) {
            x[i][0] += timeStepSize * v[i][0];
            x[i][1] += timeStepSize * v[i][1];
            x[i][2] += timeStepSize * v[i][2];
          }

          #pragma omp parallel for reduction(max: maxV)
          for (int i = 0; i < NumberOfBodies; i++) {
            v[i][0] += timeStepSize * force0[i] / mass[i];
            v[i][1] += timeStepSize * force1[i] / mass[i];
            v[i][2] += timeStepSize * force2[i] / mass[i];
            double speed = std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] );
            maxV = std::max(maxV, speed);
          }

          t += timeStepSize;

          delete[] force0;
          delete[] force1;
          delete[] force2;
        }
        double force_calculation(int i, int j, int direction) {
                    const double distance = sqrt(
                        (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
                        (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
                        (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])
                    );

                    const double distance3 = distance * distance * distance;

                    #pragma omp critical
                    minDx = std::min(minDx, distance);

                    return (x[i][direction] - x[j][direction]) * mass[i] * mass[j] / distance3;
                }
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
