#include "NBodySimulation.h"

// Need to test this

class NBodySimulationVectorised : public NBodySimulation {
    void updateBody() {
        timeStepCounter++;
        maxV = 0.0;
        minDx = std::numeric_limits<double>::max();

        // force0 = force along x direction
        // force1 = force along y direction
        // force2 = force along z direction
        double* force0 = new double[NumberOfBodies]();
        double* force1 = new double[NumberOfBodies]();
        double* force2 = new double[NumberOfBodies]();

        if (NumberOfBodies == 1) minDx = 0;  // No distances to calculate
    
        #pragma omp simd reduction(min:minDx)
        for (int i = 0; i < NumberOfBodies; i++) {
            for (int j = i + 1; j < NumberOfBodies; j++) {
                if (i != j) {
                    // x,y,z forces acting on particle i.
                    force0[i] += force_calculation(i, j, 0);
                    force1[i] += force_calculation(i, j, 1);
                    force2[i] += force_calculation(i, j, 2);
                    // x,y,z symmetric forces acting on particle j.
                    force0[j] -= force_calculation(i, j, 0);
                    force1[j] -= force_calculation(i, j, 1);
                    force2[j] -= force_calculation(i, j, 2);
                }
            }
        }

        // Position updates in a vectorised loop
        #pragma omp simd
        for (int i = 0; i < NumberOfBodies; i++) {
            x[i][0] = x[i][0] + timeStepSize * v[i][0];
            x[i][1] = x[i][1] + timeStepSize * v[i][1];
            x[i][2] = x[i][2] + timeStepSize * v[i][2];
        }

        // Velocity updates in a vectorised loop with reduction for maxV
        #pragma omp simd reduction(max:maxV)
        for (int i = 0; i < NumberOfBodies; i++) {
            v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
            v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
            v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

            maxV = std::max(maxV, std::sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]));
        }

        t += timeStepSize;

        delete[] force0;
        delete[] force1;
        delete[] force2;
    }

    double force_calculation(int i, int j, int direction) {
        #pragma omp simd
        const double distance = sqrt(
            (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
            (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
            (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])
        );

        const double distance3 = distance * distance * distance;
        minDx = std::min(minDx, distance);

        return (x[i][direction] - x[j][direction]) * mass[i] * mass[j] / distance3;
    }
};
