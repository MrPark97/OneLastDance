//
// Created by MrPark on 31.12.2019.
//
#include "task_34_08.h"

int evc_34_08(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision) {

    if(n == 1) {
        E[n-1] = A[n*n-1];
        return 0;
    }
    const int origin_n = n;
    int i = 0, j = 0, k = 0, c = 0;
    double sk = 0.0, norm_a1 = 0.0, x1 = 0.0, norm_x = 0.0;
    double scalar_y_x = 0.0, alpha = 0.0, scalar_q_x = 0.0, beta = 0.0;
    double *Q = tmp, *x = &tmp[n*n];

    // variables for speed-up modifications (shifts, exhaustions)
    double shift_k = 0.0, max_row_norm = 0.0, row_norm = 0.0;

    // calculate maximal row norm for A (||A||inf)
    for(i = 0; i < n; i++) {
        row_norm = 0.0;
        for(j = 0; j < n; j++)
            row_norm += fabs(A[i*origin_n + j]);

        max_row_norm = row_norm > max_row_norm ? row_norm : max_row_norm;
    }

    max_iterations = max_iterations > 0 ? max_iterations : n*100000;

    for(c = 0; c < max_iterations && n > 2; c++) {
        shift_k = A[origin_n*(n-1) + (n-1)];

        // A - shift * I
        for(i = 0; i < n; i++) {
            A[i*origin_n + i] -= shift_k;
            A[i*origin_n + i] = fabs(A[i*origin_n + i]) > precision ? A[i*origin_n + i] : 0;
        }

        // Q initialization
        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++)
                Q[i*origin_n + j] = 0.0;
            Q[i*origin_n + i] = 1.0;
        }

        // QR - decomposition (Householder reflection)
        for(i = 0; i < (n-1); i++) {
            sk = A[(i+1)*origin_n + i] * A[(i+1)*origin_n + i];

            if(sk <= precision)
                continue;

            norm_a1 = sqrt(A[i*origin_n + i] * A[i*origin_n + i] + sk);

            x1 = A[i*origin_n + i] - norm_a1;
            norm_x = sqrt(x1*x1 + sk);

            x[i] = x1 / norm_x;
            x[i+1] = A[(i+1)*origin_n + i] / norm_x;

            for(k = 0; k < n; k++) {
                scalar_y_x = 0.0, scalar_q_x = 0.0;
                for(j = i; j < i+2; j++) {
                    scalar_y_x += A[j*origin_n + k] * x[j];
                    scalar_q_x += Q[k*origin_n + j] * x[j];
                }

                alpha = 2.0 * scalar_y_x, beta = 2.0 * scalar_q_x;

                for(j = i; j < i+2; j++) {
                    A[j*origin_n + k] = A[j*origin_n + k] - alpha * x[j];

                    A[j*origin_n + k] = fabs(A[j*origin_n + k]) > precision ? A[j*origin_n + k] : 0.0;

                    Q[k*origin_n + j] = Q[k*origin_n + j] - beta * x[j];

                    Q[k*origin_n + j] = fabs(Q[k*origin_n + j]) > precision ? Q[k*origin_n + j] : 0.0;
                }
            }
        }

        for(i = 0; i < n; i++) {
            for(j = 0; j < n; j++)
                x[j] = A[i*origin_n + j];

            for(j = 0; j < n; j++) {
                A[i*origin_n + j] = 0.0;
                for(k = 0; k < n; k++) {
                    A[i * origin_n + j] += x[k] * Q[k*origin_n + j];
                }

                A[i*origin_n + j] = fabs(A[i*origin_n + j]) > precision ? A[i*origin_n + j] : 0.0;
            }
        }

        // A + shift * I
        for(i = 0; i < n; i++) {
            A[i*origin_n + i] += shift_k;
            A[i*origin_n + i] = fabs(A[i*origin_n + i]) > precision ? A[i*origin_n + i] : 0;
        }

        // calculate λn with epsilon
        if(fabs(A[origin_n*(n-1)+(n-2)]) <= epsilon * max_row_norm || fabs(A[origin_n*(n-1)+(n-1)]-shift_k) <= epsilon) {
            E[n-1] = A[origin_n*(n-1)+(n-1)];
            n--;
        }
    }



    if(n == 2) {
        //        |a11-λ   a12|
        //  solve |           | = 0 equation
        //        |a21   a22-λ|
        double p = - (A[origin_n*0 + 0] + A[origin_n * 1 + 1]), q = A[origin_n * 0 + 0]*A[origin_n * 1 + 1] - A[origin_n * 0 + 1]*A[origin_n * 1 + 0];
        double D = p*p - 4*q;

        if(D < -precision)
            return -1;

        if(fabs(D) > precision) {
            E[0] = (-p + sqrt(D)) / 2.0;
            E[0] = fabs(E[0]) > precision ? E[0] : 0.0;

            E[1] = (-p - sqrt(D)) / 2.0;
            E[1] = fabs(E[1]) > precision ? E[1] : 0.0;
        } else {
            E[0] = -p / 2.0;
            E[0] = fabs(E[0]) > precision ? E[0] : 0.0;
            E[1] = E[0];
        }
    } else {
        for(i = 0; i < origin_n; i++)
            E[i] = A[i*origin_n + i];
        return 1;
    }

    return 0;
}

int evc_memsize_34_08(int n) {
    return (n*n + n) * sizeof(double);
}