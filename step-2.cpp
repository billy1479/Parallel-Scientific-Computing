#include <iomanip>
#include <immintrin.h>
#include <omp.h>
#include <vector>

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

class NBodySimulationParallelised : public NBodySimulation {
    public:
        NBodySimulationParallelised() {
            omp_set_num_threads(omp_get_max_threads()); // Set number of threads to maximum available
            std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;
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

            int i = 0;
            // reduction(+:force0[:NumberOfBodies],force1[:NumberOfBodies],force2[:NumberOfBodies])
            #pragma omp parallel for schedule(dynamic) 
            for (i = 0; i<NumberOfBodies; i++) {
                #pragma omp parallel for
                for (int j=i+1; j<NumberOfBodies; j++) {
                    if(i!=j){
                        double f0 = force_calculation(i,j,0);
                        double f1 = force_calculation(i,j,1);
                        double f2 = force_calculation(i,j,2);
                        
                        // Tried atomic pragma but it was slower
                        force0[i] += f0;
                        force1[i] += f1;
                        force2[i] += f2;
                        force0[j] -= f0;
                        force1[j] -= f1;
                        force2[j] -= f2;
                    }
                }
            }

            #pragma omp parallel for shared(timeStepSize) private(i)
            for (i = 0; i < NumberOfBodies; i++){
                x[i][0] = x[i][0] + timeStepSize * v[i][0];
                x[i][1] = x[i][1] + timeStepSize * v[i][1];
                x[i][2] = x[i][2] + timeStepSize * v[i][2];
            }

            #pragma omp parallel for reduction(max:maxV) shared(timeStepSize) private(i)
            for (i = 0; i < NumberOfBodies; i++){
                v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
                v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
                v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

                maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
            }

            t += timeStepSize;

            // delete[] force0;
            // delete[] force1;
            // delete[] force2;
        }

        double force_calculation (int i, int j, int direction){
            // #pragma omp simd reduction(min:minDx)
            // Euclidean distance
            const double distance = sqrt(
                                        (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                                        (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                                        (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                                        );


            const double distance3 = distance * distance * distance;

            // #pragma omp critical(minDx) // This harms performance -> mention in report?
            // #pragma omp atomic write // This harms performance, but less than critical -> mention in report?
            minDx = std::min( minDx,distance );

            return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / distance3;
        }
    private:
        int i; // Declared as private variable to use it in shared pragmas
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
