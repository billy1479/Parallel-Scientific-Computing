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

            #pragma omp simd collapse(2) reduction(+: force0[:NumberOfBodies], force1[:NumberOfBodies], force2[:NumberOfBodies]) // HAVE OPTIMISED BY USING SIMD TO vectorise THE LOOP
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

            #pragma omp simd
            for (int i=0; i < NumberOfBodies; i++){
                v[i][0] += v[i][0] + timeStepSize * force0[i] / mass[i];
                v[i][1] += v[i][1] + timeStepSize * force1[i] / mass[i];
                v[i][2] += v[i][2] + timeStepSize * force2[i] / mass[i];

                maxV = simd_max(maxV, simd_sqrt( v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2] ));
            }

            t += timeStepSize;

            delete[] force0;
            delete[] force1;
            delete[] force2;
        }

        #pragma omp declare simd
        double force_calculation (int i, int j, int direction){
            const double distance = sqrt(
                               (x[j][0]-x[i][0]) * (x[j][0]-x[i][0]) +
                               (x[j][1]-x[i][1]) * (x[j][1]-x[i][1]) +
                               (x[j][2]-x[i][2]) * (x[j][2]-x[i][2])
                               );

            minDx = simd_min( minDx,distance);

            return (x[i][direction]-x[j][direction]) * mass[i]*mass[j] / (distance * distance * distance);
        }
};
