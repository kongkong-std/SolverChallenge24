/*
 * This is the reference direct method solver code of SolverChallenge24
 * author: Li Zhao and Haibing Sun
 */

#include "mysolver.h"

// read matrix file
#include "file_process.h"

// utlise function file
#include "utlise.h"
#include "utlise_long.h"

int main(int argc, char **argv)
{
    int m, n, nnzA, isSymmetricA;
    int *row_ptr = NULL;            // the csr row pointer array of matrix A
    int *col_idx = NULL;            // the csr column index array of matrix A
    double *val = NULL;             // the csr value array of matrix A (real number)
    double *val_im = NULL;          // the csr value array of matrix A (imaginary number)
    double *x = NULL, *x_im = NULL; // solution vector x, (x: real number, x_im: imaginary number)
    double *b = NULL, *b_im = NULL; // right-hand side vector b, (b: real number, b_im: imaginary number)
    double tt, time;

    char *filename_matrix;    // the filename of matrix A
    char *filename_b;         // the filename of right-hand side vector b
    int read_matrix_base = 1; // 0-base or 1-base, default 1-base
    int type = 0;             // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    int test_frequency = 1;   // run code frequency
    int sys_type = 0;         // type of algebraic systems, 0: real, 1: complex; default 0

    /* ========================================== */
    // Step 0: Read command line argument
    /* ========================================== */
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-file_mat", argv[index]))
        {
            filename_matrix = argv[index + 1];
        }
        if (strstr("-file_rhs", argv[index]))
        {
            filename_b = argv[index + 1];
        }
        if (strstr("-type", argv[index]))
        {
            type = atoi(argv[index + 1]);
        }
        if (strstr("-sys_type", argv[index]))
        {
            sys_type = atoi(argv[index + 1]);
        }
        if (strstr("-test_frequency", argv[index]))
        {
            test_frequency = atoi(argv[index + 1]);
        }
    }

    fprintf(stdout, "matrix name :      %s\nvectorb name :     %s\nread function :    base-%d\ntype :             %d\nsys_type :         %d\n",
            filename_matrix, filename_b, read_matrix_base, type, sys_type);

    /* ========================================== */
    // Step 1: Load matrix and rhs from mtx files
    /* ========================================== */
    if (sys_type == 0) // real system
    {
        printf("\n>>>> begin to process linear system file...\n");
        RealCOO2CSRMatrixFileProcess(filename_matrix, &m, &n, &nnzA, &row_ptr, &col_idx, &val);

        if (m != n)
        {
            fprintf(stdout, "Invalid matrix size.\n");
            exit(EXIT_FAILURE);
        }

        x = (double *)malloc(sizeof(double) * n);

        // load right-hand side vector b
        RealRHSFileProcess(filename_b, &b);

        // initial vector x
        memset(x, 0.0, sizeof(double) * n);
        printf(">>>> linear system file has been processed!!!\n\n");

        printf(">>>> begin to call kml-dss solver...\n");
        // direct solver sample
        MySolver mysolver;
        KMLRealSolverInitialize(&mysolver, x, n, row_ptr, col_idx, val, b);

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLRealSolverAnalyze(&mysolver);
                KMLRealSolverFactor(&mysolver);
                KMLRealSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            KMLRealSolverAnalyze(&mysolver);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLRealSolverFactor(&mysolver);
                KMLRealSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time
            KMLRealSolverAnalyze(&mysolver);
            KMLRealSolverFactor(&mysolver);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLRealSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

        KMLRealSolverQuery(&mysolver);

        // free kml-dss handle
        KMLRealSolverClean(&mysolver);
        printf(">>>> kml-dss solver has done!!!\n\n");

        // check time
        if (type == 0)
        {
            fprintf(stdout, "CHECK end to end time :         %12.6lf ms\n", time);
        }
        else if (type == 1)
        {
            fprintf(stdout, "CHECK solver + solve time :     %12.6lf ms\n", time);
        }
        else if (type == 2)
        {
            fprintf(stdout, "CHECK solve time :              %12.6lf ms\n", time);
        }

        // check the memory
        mem_usage();

        // 1) using double precision
        check_correctness(n, row_ptr, col_idx, val, x, b);
        // 2) using long double precision
        check_correctness_ld_d2ld(n, row_ptr, col_idx, val, x, b);

        free(row_ptr);
        free(col_idx);
        free(val);
        free(x);
        free(b);
    }

    if (sys_type == 1)
    { // complex system
        printf("\n>>>> begin to process linear system file...\n");
        ComplexCOO2CSRMatrixFileProcess(filename_matrix, &m, &n, &nnzA, &row_ptr, &col_idx, &val, &val_im);

        if (m != n)
        {
            fprintf(stdout, "Invalid matrix size.\n");
            return 0;
        }

        x = (double *)malloc(sizeof(double) * n);
        x_im = (double *)malloc(sizeof(double) * n);
        b = (double *)malloc(sizeof(double) * n);
        b_im = (double *)malloc(sizeof(double) * n);

        ComplexRHSFileProcess(filename_b, &b, &b_im);

        // initial vector x
        memset(x, 0.0, sizeof(double) * n);
        memset(x_im, 0.0, sizeof(double) * n);
        printf(">>>> linear system file has been processed!!!\n\n");

        printf(">>>> begin to call kml-dss solver...\n");
        MySolverComplex mysolver;
        kml_complex_double *kml_dss_solver_x = NULL, *kml_dss_solver_val = NULL, *kml_dss_solver_b = NULL;
        if ((kml_dss_solver_x = (kml_complex_double *)malloc(n * sizeof(kml_complex_double))) == NULL ||
            (kml_dss_solver_b = (kml_complex_double *)malloc(n * sizeof(kml_complex_double))) == NULL ||
            (kml_dss_solver_val = (kml_complex_double *)malloc(nnzA * sizeof(kml_complex_double))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - kml_complex_double\n");
            exit(EXIT_FAILURE);
        }

        for (int index = 0; index < n; ++index)
        {
            kml_dss_solver_x[index] = x[index] + x_im[index] * _Complex_I;
            kml_dss_solver_b[index] = b[index] + b_im[index] * _Complex_I;
        }
        for (int index = 0; index < nnzA; ++index)
        {
            kml_dss_solver_val[index] = val[index] + val_im[index] * _Complex_I;
        }

        KMLComplexSolverInitialize(&mysolver, kml_dss_solver_x, n,
                                   row_ptr, col_idx, kml_dss_solver_val, kml_dss_solver_b);

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLComplexSolverAnalyze(&mysolver);
                KMLComplexSolverFactor(&mysolver);
                KMLComplexSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            KMLComplexSolverAnalyze(&mysolver);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLComplexSolverFactor(&mysolver);
                KMLComplexSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time

            KMLComplexSolverAnalyze(&mysolver);
            KMLComplexSolverFactor(&mysolver);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                KMLComplexSolverSolve(&mysolver);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

        KMLComplexSolverQuery(&mysolver);

        for (int index = 0; index < n; ++index)
        {
            x[index] = creal(kml_dss_solver_x[index]);
            x_im[index] = cimag(kml_dss_solver_x[index]);
        }

        free(kml_dss_solver_val);
        free(kml_dss_solver_b);
        free(kml_dss_solver_x);

        // free kml-dss handle
        KMLComplexSolverClean(&mysolver);
        printf(">>>> kml-dss solver has done!!!\n\n");

        // check time
        if (type == 0)
        {
            fprintf(stdout, "CHECK end to end time :         %12.6lf ms\n", time);
        }
        else if (type == 1)
        {
            fprintf(stdout, "CHECK solver + solve time :     %12.6lf ms\n", time);
        }
        else if (type == 2)
        {
            fprintf(stdout, "CHECK solve time :              %12.6lf ms\n", time);
        }

        // check the memory
        mem_usage();

        // 1) using double precision
        check_correctness_complex(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);
        // 2) using long double precision
        check_correctness_complex_ld_d2ld(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);

        free(row_ptr);
        free(col_idx);
        free(val);
        free(x);
        free(b);
        free(val_im);
        free(x_im);
        free(b_im);
    }

    fprintf(stdout, "------------------------------------------\n");

    /* ========================================== */
    // Step 3: Check time, memory and correctness
    /* ========================================== */

#if 0
    // store x to a file
    char *answer_x = "answer_x.rhs";
    if (sys_type == 0) // real rhs
        store_x(n, x, answer_x);
    else // complex rhs
        store_x_complex(n, x, x_im, answer_x);
#endif // store solution file

    return 0;
}
