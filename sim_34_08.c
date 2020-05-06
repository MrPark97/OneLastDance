//
// Created by MrPark on 31.12.2019.
//
#include "task_34_08.h"

int sim_34_08(int n, double* A, double* tmp, double precision) {
    int i=0, j=0, k=0;
    double cosine = 0.0, sine = 0.0, x = 0.0, y = 0.0, xi = 0.0, xj = 0.0;
    for(j = 0; j < n; j++) {
        for(i = j+2; i < n; i++) {
            x = A[(j+1)*n + j];
            y = A[i*n+j];

            if(fabs(x) <= precision && fabs(y) <= precision) {
                continue;
            }

            cosine = x / sqrt(x*x + y*y);
            sine = - y / sqrt(x*x + y*y);

            // T * A
            for(k = 0; k < n; k++) {
                xi = A[(j+1)*n + k];
                xj = A[i*n + k];

                A[(j+1)*n + k] = xi * cosine - xj * sine;
                A[(j+1)*n + k] = fabs(A[(j+1)*n + k]) > precision ? A[(j+1)*n + k] : 0.0;

                A[i*n + k] = xi * sine + xj * cosine;
                A[i*n + k] = fabs(A[i*n + k]) > precision ? A[i*n + k] : 0.0;
            }

            // T * A * Tt
            for(k = 0; k  < n; k++) {

                xi = A[k*n + (j+1)];
                xj = A[k*n + i];

                A[k*n + (j+1)] = xi * cosine - xj * sine;
                A[k*n + (j+1)] = fabs(A[k*n + (j+1)]) > precision ? A[k*n + (j+1)] : 0.0;

                A[k*n + i] = xi * sine + xj * cosine;
                A[k*n + i] = fabs(A[k*n + i]) > precision ? A[k*n + i] : 0.0;
            }
        }
    }

    return 0;
}



int sim_memsize_34_08(int n) {
    return 0;
}

