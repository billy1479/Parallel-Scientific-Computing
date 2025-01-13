#include "NBodySimulation.h"
#include <iostream>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cstdlib>

class NBodySimulationVectorised : public NBodySimulation {
    public:
        void updateBody () {
            timeStepCounter++;
            maxV   = 0.0;
            minDx  = std::numeric_limits<double>::max();

            double* force0 = new double[NumberOfBodies]();
            double* force1 = new double[NumberOfBodies]();
            double* force2 = new double[NumberOfBodies]();

            if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate

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
                x[i][0] += timeStepSize * v[i][0];
                x[i][1] += timeStepSize * v[i][1];
                x[i][2] += timeStepSize * v[i][2];
            }

            for (int i=0; i < NumberOfBodies; i++){
                v[i][0] += v[i][0] + timeStepSize * force0[i] / mass[i];
                v[i][1] += v[i][1] + timeStepSize * force1[i] / mass[i];
                v[i][2] += v[i][2] + timeStepSize * force2[i] / mass[i];

                maxV = std::max(maxV, std::sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
            }

            t += timeStepSize;

            delete[] force0;
            delete[] force1;
            delete[] force2;
        }

        #pragma omp declare simd
        double force_calculation (int i, int j, int direction){// HAVE OPTIMISED BY USING VARIABLES TO STORE THE FORCE CALCULATION                            
            double temp_1  = x[j][0]-x[i][0];
            double temp_2  = x[j][1]-x[i][1];
            double temp_3  = x[j][2]-x[i][2];

            double distance = sqrt( temp_1 * temp_1 + temp_2 * temp_2 + temp_3 * temp_3 );

            minDx = std::min( minDx,distance);

            return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / (distance * distance * distance);
        }
};
