//
// Created by MrPark on 31.12.2019.
//

#ifndef ONELASTDANCE_TASK_34_08_H
#define ONELASTDANCE_TASK_34_08_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAX_STRING_LENGTH 256

int sim_34_08(int n, double* A, double* tmp, double precision);

int sim_memsize_34_08(int n);

int evc_34_08(int n, int max_iterations, double epsilon, double* A, double* E, double* tmp, double precision);

int evc_memsize_34_08(int n);

#endif //ONELASTDANCE_TASK_34_08_H
