#include "task_34_08.h"

char error_msgs[][MAX_STRING_LENGTH] = {
        "",
        "",
        "using -h, -? options error", // begins from 2 because of task format
        "undefined option error",
        "wrong interface parameters error",
        "cannot open output file error",
        "cannot open input file error",
        "reading data error"
};

// comparator for qsort
int cmp(const void *x, const void *y)
{
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}


int main(int argc, char* argv[])
{
    clock_t start_time = clock();
    double precision = 1e-14, epsilon = 1e-10;
    char *input_filename = "34_08_in.txt", *output_filename = "34_08_out.txt";

    int i = 0, max_iter = 0;

    char FLAG_D = 0, FLAG_E = 0, FLAG_P = 0, FLAG_T = 0; // flags (cmd options)
    char FLAG_F = 0; // custom filenames

    for (i = 1; i < argc; ++i)
    {
        if (argv[i][0] == '-' && argv[i][2] == '\0')
            switch (argv[i][1])
            {
                case 'd':
                    FLAG_D = 1;
                    break;
                case 'e':
                    FLAG_E = 1;
                    break;
                case 'p':
                    FLAG_P = 1;
                    break;
                case 't':
                    FLAG_T = 1;
                    break;
                case 'h':;
                case '?':
                    printf("Usage: lss [input_file_name] [output_file_name] [options]\n"
                           "Where options include:\n"
                           "-d    print debug messages [default OFF]\n"
                           "-e    print errors [default OFF]\n"
                           "-p    print matrix [default OFF]\n"
                           "-t    print execution t [default OFF]\n"
                           "-prec=<num>       precision [default - 1e-14]\n"
                           "-eps=<num>        epsilon [default - 1e-10]"
                           "-max_iter=<num>   limit number of iterations [default - 0, i.e. not limit]"
                           "-h, -?     print this and exit\n");
                    return 2;
                default:
                    return 3;
            }
        else if(argv[i][0] == '-') {
            if(argv[i][1] == 'p' && argv[i][2] == 'r' && argv[i][3] == 'e' && argv[i][4] == 'c' && argv[i][5] == '=') {
                precision = atof(&argv[i][6]);
            } else if(argv[i][1] == 'e' && argv[i][2] == 'p' && argv[i][3] == 's' && argv[i][4] == '=') {
                epsilon = atof(&argv[i][5]);
            } else if(argv[i][1] == 'm' && argv[i][2] == 'a' && argv[i][3] == 'x' && argv[i][4] == '_' && argv[i][5] == 'i' && argv[i][6] == 't' && argv[i][7] == 'e' && argv[i][8] == 'r' && argv[i][9] == '=') {
                max_iter = atoi(&argv[i][10]);
            } else {
                return 3;
            }
        }
        else
        {
            if (i == 1)
                input_filename = argv[i], FLAG_F = 1;
            else if (i == 2 && FLAG_F)
                output_filename = argv[i];
            else
            {
                return 4;
            }
        }
    }

    FILE *input_file, *output_file;

    if ((output_file = fopen(output_filename, "w")) == NULL)
    {
        if(FLAG_E)
            fprintf(stderr, "%s\n", error_msgs[5]);
        return 5;
    }

    if ((input_file = fopen(input_filename, "r")) == NULL)
    {
        if(FLAG_E)
            fprintf(stderr, "%s\n", error_msgs[6]);
        fclose(output_file);
        return 6;
    }

    int n = 0;

    if (fscanf(input_file, "%d", &n) != 1)
    {
        if(FLAG_E)
            fprintf(stderr, "%s\n", error_msgs[7]);
        fclose(output_file);
        fclose(input_file);
        return 7;
    }

    double *A = (double *) malloc(n * n * sizeof(double));

    int j = 0;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (fscanf(input_file, "%lf", &A[i * n + j]) != 1)
            {
                if(FLAG_E)
                    fprintf(stderr, "%s\n", error_msgs[7]);
                fclose(output_file);
                fclose(input_file);
                return 7;
            }

    if (FLAG_P)
    {
        printf("%d\n", n);

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                printf("%1.9lf ", A[i * n + j]);
            printf("\n");
        }

        printf("\n");
    }

    double *tmp = (double *) malloc(sim_memsize_34_08(n));

    int status = sim_34_08(n, A, tmp, precision);

    if(status == -1) {
        fprintf(output_file, "%d\n", 0);
        return -1;
    }

    tmp = (double *) malloc(evc_memsize_34_08(n));
    double *E = (double *) malloc(n * sizeof(double));

    status = evc_34_08(n, max_iter, epsilon, A, E, tmp, precision);

    if(status == -1) {
        fprintf(output_file, "%d\n", 0);
        return -1;
    }

    clock_t end_time = clock();
    double execution_time = (double)((end_time - start_time)) / CLOCKS_PER_SEC;

    if (FLAG_T)
        printf("Execution time: %1.9lf \n", execution_time);

    if(status == 0) {
        qsort(E, n, sizeof(double), cmp);
    }

    if(status == 0 || status == 1)
    {
        fprintf(output_file, "%d\n", n);

        for (i = 0; i < n; ++i)
            fprintf(output_file, "%1.9lf\n", E[i]);
    }

    fclose(input_file);
    fclose(output_file);
    free(A);
    free(E);
    free(tmp);
    return status;
}