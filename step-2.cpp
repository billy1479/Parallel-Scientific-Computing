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

class NBodySimulationParallelised : public NBodySimulation {
    public:
        NBodySimulationParallelised() {
            omp_set_num_threads(omp_get_max_threads()); // Set number of threads to maximum available
            omp_set_nested(0); // Disable nested parallelism -> reduces overhead
            std::cout << "Number of threads: " << omp_get_max_threads() << std::endl;   
        }

        void setUp (int argc, char** argv) {
            checkInput(argc, argv);

            NumberOfBodies = (argc-4) / 7;

            x0 = new double[NumberOfBodies];  // x direction positions
            x1 = new double[NumberOfBodies];  // y direction positions
            x2 = new double[NumberOfBodies];  // z direction positions
            
            v0 = new double[NumberOfBodies];  // x direction velocities
            v1 = new double[NumberOfBodies];  // y direction velocities
            v2 = new double[NumberOfBodies];  // z direction velocities
            
            mass = new double[NumberOfBodies];

            int readArgument = 1;

            tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
            tFinal       = std::stof(argv[readArgument]); readArgument++;
            timeStepSize = std::stof(argv[readArgument]); readArgument++;

            for (int i=0; i<NumberOfBodies; i++) {
                // Modified array declarations

                // Reading positions
                x0[i] = std::stof(argv[readArgument]); readArgument++;
                x1[i] = std::stof(argv[readArgument]); readArgument++;
                x2[i] = std::stof(argv[readArgument]); readArgument++;

                // Reading velocities
                v0[i] = std::stof(argv[readArgument]); readArgument++;
                v1[i] = std::stof(argv[readArgument]); readArgument++;
                v2[i] = std::stof(argv[readArgument]); readArgument++;

                mass[i] = std::stof(argv[readArgument]); readArgument++;
                
                if (mass[i]<=0.0 ) {
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

            double* force0 = new double[NumberOfBodies]();
            double* force1 = new double[NumberOfBodies]();
            double* force2 = new double[NumberOfBodies]();

            if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate
    
            #pragma omp parallel for \
                reduction(min:minDx) \
                reduction(+:force0[:NumberOfBodies]) \
                reduction(+:force1[:NumberOfBodies]) \
                reduction(+:force2[:NumberOfBodies])
            for (int i = 0; i < NumberOfBodies - 1; i++) {
                for (int j = i + 1; j < NumberOfBodies; j++) {
                    const double dx = x0[j] - x0[i];
                    const double dy = x1[j] - x1[i];
                    const double dz = x2[j] - x2[i];

                    const double distance = sqrt(dx*dx + dy*dy + dz*dz);
                    minDx = std::min(minDx, distance);

                    const double factor = mass[i] * mass[j] / (distance * distance * distance);
                    
                    const double fx = dx * factor;
                    const double fy = dy * factor;
                    const double fz = dz * factor;

                    // Direct updates to force arrays - no local arrays needed
                    force0[i] -= fx;
                    force1[i] -= fy;
                    force2[i] -= fz;
                    force0[j] += fx;
                    force1[j] += fy;
                    force2[j] += fz;
                }
            }
            
            #pragma omp parallel for simd schedule(static) shared(timeStepSize)
            for (i = 0; i < NumberOfBodies; i++){
                x0[i] = x0[i] + timeStepSize * v0[i];
                x1[i] = x1[i] + timeStepSize * v1[i];
                x2[i] = x2[i] + timeStepSize * v2[i];

                v0[i] = v0[i] + timeStepSize * force0[i] / mass[i];
                v1[i] = v1[i] + timeStepSize * force1[i] / mass[i];
                v2[i] = v2[i] + timeStepSize * force2[i] / mass[i];

                maxV = std::max(maxV, std::sqrt(v0[i]*v0[i] + v1[i]*v1[i] + v2[i]*v2[i]));
            }

            t += timeStepSize;

            delete[] force0;
            delete[] force1;
            delete[] force2;
        }

        bool hasReachedEnd () {
            return t > tFinal;
            }

            void takeSnapshot () {
            if (t >= tPlot) {
                printParaviewSnapshot();
                printSnapshotSummary();
                tPlot += tPlotDelta;
            }
            }

            void openParaviewVideoFile () {
            videoFile.open("paraview-output/result.pvd");
            videoFile << "<?xml version=\"1.0\"?>" << std::endl
                        << "<VTKFile type=\"Collection\""
                " version=\"0.1\""
                " byte_order=\"LittleEndian\""
                " compressor=\"vtkZLibDataCompressor\">" << std::endl
                        << "<Collection>";
            }

            void closeParaviewVideoFile () {
            videoFile << "</Collection>"
                        << "</VTKFile>" << std::endl;
            videoFile.close();
            }

            void printParaviewSnapshot () {
            static int counter = -1;
            counter++;
            std::stringstream filename, filename_nofolder;
            filename << "paraview-output/result-" << counter <<  ".vtp";
            filename_nofolder << "result-" << counter <<  ".vtp";
            std::ofstream out( filename.str().c_str() );
            out << "<VTKFile type=\"PolyData\" >" << std::endl
                << "<PolyData>" << std::endl
                << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
                << "  <Points>" << std::endl
                << "   <DataArray type=\"Float64\""
                " NumberOfComponents=\"3\""
                " format=\"ascii\">";

            for (int i=0; i<NumberOfBodies; i++) {
                out << x0[i]
                    << " "
                    << x1[i]
                    << " "
                    << x2[i]
                    << " ";
            }

            out << "   </DataArray>" << std::endl
                << "  </Points>" << std::endl
                << " </Piece>" << std::endl
                << "</PolyData>" << std::endl
                << "</VTKFile>"  << std::endl;

            out.close();

            videoFile << "<DataSet timestep=\"" << counter
                        << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
                        << "\"/>" << std::endl;
            }

            void printSnapshotSummary () {
            std::cout << "plot next snapshot"
                        << ",\t time step=" << timeStepCounter
                        << ",\t t="         << t
                        << ",\t dt="        << timeStepSize
                        << ",\t v_max="     << maxV
                        << ",\t dx_min="    << minDx
                        << std::endl;
            }

            void printSummary () {
            std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
            std::cout << "Position of first remaining object: "
                        << x0[0] << ", " << x1[0] << ", " << x2[0] << std::endl;
            }
                private:
                    int i; // Declared as private variable to use it in shared
                    double *x0, *x1, *x2, *v0, *v1, *v2, *mass;
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
