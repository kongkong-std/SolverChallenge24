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
    int m, n, nnzA;
    int *row_ptr = NULL;            // the csr row pointer array of matrix A
    int *col_idx = NULL;            // the csr column index array of matrix A
    double *val = NULL;             // the csr value array of matrix A (real number)
    double *val_im = NULL;          // the csr value array of matrix A (imaginary number)
    double *x = NULL, *x_im = NULL; // solution vector x, (x: real number, x_im: imaginary number)
    double *b = NULL, *b_im = NULL; // right-hand side vector b, (b: real number, b_im: imaginary number)
    double tt, time;
    double tt_analyze, time_analyze; // time of analyze
    double tt_factor, time_factor;   // time of factor
    double tt_solve, time_solve;     // time of solve

    char *filename_matrix;          // the filename of matrix A
    char *filename_b;               // the filename of right-hand side vector b
    int read_matrix_base = 1;       // 0-base or 1-base, default 1-base
    int type = 0;                   // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    int test_frequency = 1;         // run code frequency
    int sys_type = 0;               // type of algebraic systems, 0: real, 1: complex; default 0
    int IR_times = 10;              // max IR time
    double sys_rtol = 1e-8;         // relative residual tolerance
    int n_thread = 32;              // number of threads
    int n_thread_rdr = 8;           // number of threads of fill-in reduction
    double factor_threshold = 1e-8; // factor threshold
    char *matrix_type = NULL;       // matrix type KmlSolverMatrixType
    char *refine_method = NULL;     // solver refine method
    int refine_maxit = 1;           // refine max iter
    double refine_tol = sys_rtol;   // refine tolerance
    KmlSolverMatrixType kml_matrix_type = KMLSS_MATRIX_GEN;
    KmlDssRefineMethod kml_refine_method = KMLDSS_REFINE_OFF;

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
        if (strstr("-ir_times", argv[index]))
        {
            IR_times = atoi(argv[index + 1]);
        }
        if (strstr("-sys_rtol", argv[index]))
        {
            sys_rtol = atof(argv[index + 1]);
        }
        if ((strstr("-n_thread", argv[index])))
        {
            n_thread = atoi(argv[index + 1]);
        }
        if ((strstr("-nt_rdr", argv[index])))
        {
            n_thread_rdr = atoi(argv[index + 1]);
        }
        if (strstr("-factor_threshold", argv[index]))
        {
            factor_threshold = atof(argv[index + 1]);
        }
        if (strstr("-matrix_type", argv[index]))
        {
            matrix_type = argv[index + 1];
        }
        if (strstr("-refine_method", argv[index]))
        {
            refine_method = argv[index + 1];
        }
        if (strstr("-refine_maxit", argv[index]))
        {
            refine_maxit = atoi(argv[index + 1]);
        }
        if (strstr("-refine_tol", argv[index]))
        {
            refine_tol = atof(argv[index + 1]);
        }
    }

    fprintf(stdout, "matrix name :      %s\nvectorb name :     %s\nread function :    base-%d\ntype :             %d\nsys_type :         %d\n",
            filename_matrix, filename_b, read_matrix_base, type, sys_type);

    kml_matrix_type = ParseMatrixType(matrix_type);
    kml_refine_method = ParseRefineMethod(refine_method);

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
        KMLRealSolverMatrixCreate(&mysolver, kml_matrix_type, n, row_ptr, col_idx, val);
        KMLRealSolverRHSCreate(&mysolver, n, b);
        KMLRealSolverSOLCreate(&mysolver, n, x);
        KMLRealSolverInitialize(&mysolver, n_thread);

        // residual and error definition
        /*
         * 1. solver_r = b - Ax
         * 2. A solver_e = solver_r
         * 3. x = x + solver_e
         */
        double *solver_r = NULL, *solver_e = NULL;
        if ((solver_r = (double *)malloc(sizeof(double) * n)) == NULL ||
            (solver_e = (double *)malloc(sizeof(double) * n)) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - \'residual and error vector\'\n");
            exit(EXIT_FAILURE);
        }

        // initialize solver_e
        memset(solver_e, 0.0, sizeof(double) * n);

        double solver_r_l2_norm = 0., solver_b_l2_norm = 0.;
        solver_b_l2_norm = vec2norm(b, n);
        solver_r_l2_norm = solver_b_l2_norm; // initialize || solver_r ||_2
        int ir_times = 0;

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_analyze = GetCurrentTime();
                KMLRealSolverAnalyze(&mysolver, n_thread_rdr);
                time_analyze = GetCurrentTime() - tt_analyze;
                printf("---- time of analyze: %12.6lf ms\n", time_analyze);

                tt_factor = GetCurrentTime();
                KMLRealSolverFactor(&mysolver, factor_threshold);
                time_factor = GetCurrentTime() - tt_factor;
                printf("---- time of factor: %12.6lf ms\n", time_factor);

                tt_solve = GetCurrentTime();
                KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                spmv(n, row_ptr, col_idx, val, x, solver_r); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r[index] = b[index] - solver_r[index]; // solver_r = b - Ax
                }
                solver_r_l2_norm = vec2norm(solver_r, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    KMLRealSolverRHSCreate(&mysolver, n, solver_r);
                    KMLRealSolverSOLCreate(&mysolver, n, solver_e);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, solver_e);
                    KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        x[index] += solver_e[index];
                    }

                    spmv(n, row_ptr, col_idx, val, x, solver_r);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r[index] = b[index] - solver_r[index];
                    }
                    solver_r_l2_norm = vec2norm(solver_r, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLRealSolverRHSCreate(&mysolver, n, b);
                KMLRealSolverSOLCreate(&mysolver, n, x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            tt_analyze = GetCurrentTime();
            KMLRealSolverAnalyze(&mysolver, n_thread_rdr);
            time_analyze = GetCurrentTime() - tt_analyze;
            printf("---- time of analyze: %12.6lf ms\n", time_analyze);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_factor = GetCurrentTime();
                KMLRealSolverFactor(&mysolver, factor_threshold);
                time_factor = GetCurrentTime() - tt_factor;
                printf("---- time of factor: %12.6lf ms\n", time_factor);

                tt_solve = GetCurrentTime();
                KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                spmv(n, row_ptr, col_idx, val, x, solver_r); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r[index] = b[index] - solver_r[index]; // solver_r = b - Ax
                }
                solver_r_l2_norm = vec2norm(solver_r, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    KMLRealSolverRHSCreate(&mysolver, n, solver_r);
                    KMLRealSolverSOLCreate(&mysolver, n, solver_e);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, solver_e);
                    KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        x[index] += solver_e[index];
                    }

                    spmv(n, row_ptr, col_idx, val, x, solver_r);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r[index] = b[index] - solver_r[index];
                    }
                    solver_r_l2_norm = vec2norm(solver_r, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLRealSolverRHSCreate(&mysolver, n, b);
                KMLRealSolverSOLCreate(&mysolver, n, x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time
            tt_analyze = GetCurrentTime();
            KMLRealSolverAnalyze(&mysolver, n_thread_rdr);
            time_analyze = GetCurrentTime() - tt_analyze;
            printf("---- time of analyze: %12.6lf ms\n", time_analyze);

            tt_factor = GetCurrentTime();
            KMLRealSolverFactor(&mysolver, factor_threshold);
            time_factor = GetCurrentTime() - tt_factor;
            printf("---- time of factor: %12.6lf ms\n", time_factor);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_solve = GetCurrentTime();
                KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                spmv(n, row_ptr, col_idx, val, x, solver_r); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r[index] = b[index] - solver_r[index]; // solver_r = b - Ax
                }
                solver_r_l2_norm = vec2norm(solver_r, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    KMLRealSolverRHSCreate(&mysolver, n, solver_r);
                    KMLRealSolverSOLCreate(&mysolver, n, solver_e);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, solver_e);
                    KMLRealSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        x[index] += solver_e[index];
                    }

                    spmv(n, row_ptr, col_idx, val, x, solver_r);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r[index] = b[index] - solver_r[index];
                    }
                    solver_r_l2_norm = vec2norm(solver_r, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLRealSolverRHSCreate(&mysolver, n, b);
                KMLRealSolverSOLCreate(&mysolver, n, x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

        // free ir
        free(solver_r);
        free(solver_e);

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
        // b = (double *)malloc(sizeof(double) * n);
        // b_im = (double *)malloc(sizeof(double) * n);

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

        KMLComplexSolverMatrixCreate(&mysolver, kml_matrix_type, n, row_ptr, col_idx, kml_dss_solver_val);
        KMLComplexSolverRHSCreate(&mysolver, n, kml_dss_solver_b);
        KMLComplexSolverSOLCreate(&mysolver, n, kml_dss_solver_x);
        KMLComplexSolverInitialize(&mysolver, n_thread);

        // residual and error definition
        /*
         * 1. solver_r = b - Ax
         * 2. A solver_e = solver_r
         * 3. x = x + solver_e
         */
        double *solver_r_re = NULL, *solver_r_im = NULL;
        double *solver_e_re = NULL, *solver_e_im = NULL;

        if ((solver_r_re = (double *)malloc(n * sizeof(double))) == NULL ||
            (solver_r_im = (double *)malloc(n * sizeof(double))) == NULL ||
            (solver_e_re = (double *)malloc(n * sizeof(double))) == NULL ||
            (solver_e_im = (double *)malloc(n * sizeof(double))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - \'residual and error vector\'\n");
            exit(EXIT_FAILURE);
        }

        // initialize solver_e
        memset(solver_e_re, 0., n * sizeof(double));
        memset(solver_e_im, 0., n * sizeof(double));

        double solver_r_l2_norm = 0., solver_b_l2_norm = 0.;
        double solver_r_l2_norm_i = 0., solver_b_l2_norm_i = 0.;
        vec2norm_complex(b, b_im, &solver_b_l2_norm, &solver_b_l2_norm_i, n);
        solver_r_l2_norm = solver_b_l2_norm;     // initialize || solver_r ||_2
        solver_r_l2_norm_i = solver_b_l2_norm_i; // initialize || solver_r ||_2
        int ir_times = 0;

        kml_complex_double *kml_dss_solver_r = NULL, *kml_dss_solver_e = NULL;
        if ((kml_dss_solver_r = (kml_complex_double *)malloc(n * sizeof(kml_complex_double))) == NULL ||
            (kml_dss_solver_e = (kml_complex_double *)malloc(n * sizeof(kml_complex_double))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - \'kmlcomplexdouble residual and error vector\'\n");
            exit(EXIT_FAILURE);
        }

        for (int index = 0; index < n; ++index)
        {
            kml_dss_solver_e[index] = solver_e_re[index] + solver_e_im[index] * _Complex_I;
        }

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_analyze = GetCurrentTime();
                KMLComplexSolverAnalyze(&mysolver, n_thread_rdr);
                time_analyze = GetCurrentTime() - tt_analyze;
                printf("---- time of analyze: %12.6lf ms\n", time_analyze);

                tt_factor = GetCurrentTime();
                KMLComplexSolverFactor(&mysolver, factor_threshold);
                time_factor = GetCurrentTime() - tt_factor;
                printf("---- time of factor: %12.6lf ms\n", time_factor);

                tt_solve = GetCurrentTime();
                KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                for (int index = 0; index < n; ++index)
                {
                    x[index] = creal(kml_dss_solver_x[index]); // updating solution
                    x_im[index] = cimag(kml_dss_solver_x[index]);
                }
                spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r_re[index] = b[index] - solver_r_re[index]; // solver_r = b - Ax
                    solver_r_im[index] = b_im[index] - solver_r_im[index];
                }
                vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_r[index] = solver_r_re[index] + solver_r_im[index] * _Complex_I;
                        kml_dss_solver_e[index] = solver_e_re[index] + solver_e_im[index] * _Complex_I;
                    }
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_e);
                    KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_x[index] += kml_dss_solver_e[index];
                        x[index] = creal(kml_dss_solver_x[index]);
                        x_im[index] = cimag(kml_dss_solver_x[index]);
                    }

                    spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r_re[index] = b[index] - solver_r_re[index];
                        solver_r_im[index] = b_im[index] - solver_r_im[index];
                    }
                    vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLComplexSolverRHSCreate(&mysolver, n, kml_dss_solver_b);
                KMLComplexSolverSOLCreate(&mysolver, n, kml_dss_solver_x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            tt_analyze = GetCurrentTime();
            KMLComplexSolverAnalyze(&mysolver, n_thread_rdr);
            time_analyze = GetCurrentTime() - tt_analyze;
            printf("---- time of analyze: %12.6lf ms\n", time_analyze);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_factor = GetCurrentTime();
                KMLComplexSolverFactor(&mysolver, factor_threshold);
                time_factor = GetCurrentTime() - tt_factor;
                printf("---- time of factor: %12.6lf ms\n", time_factor);

                tt_solve = GetCurrentTime();
                KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                for (int index = 0; index < n; ++index)
                {
                    x[index] = creal(kml_dss_solver_x[index]); // updating solution
                    x_im[index] = cimag(kml_dss_solver_x[index]);
                }
                spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r_re[index] = b[index] - solver_r_re[index]; // solver_r = b - Ax
                    solver_r_im[index] = b_im[index] - solver_r_im[index];
                }
                vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_r[index] = solver_r_re[index] + solver_r_im[index] * _Complex_I;
                        kml_dss_solver_e[index] = solver_e_re[index] + solver_e_im[index] * _Complex_I;
                    }
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_e);
                    KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_x[index] += kml_dss_solver_e[index];
                        x[index] = creal(kml_dss_solver_x[index]);
                        x_im[index] = cimag(kml_dss_solver_x[index]);
                    }

                    spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r_re[index] = b[index] - solver_r_re[index];
                        solver_r_im[index] = b_im[index] - solver_r_im[index];
                    }
                    vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLComplexSolverRHSCreate(&mysolver, n, kml_dss_solver_b);
                KMLComplexSolverSOLCreate(&mysolver, n, kml_dss_solver_x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time
            tt_analyze = GetCurrentTime();
            KMLComplexSolverAnalyze(&mysolver, n_thread_rdr);
            time_analyze = GetCurrentTime() - tt_analyze;
            printf("---- time of analyze: %12.6lf ms\n", time_analyze);

            tt_factor = GetCurrentTime();
            KMLComplexSolverFactor(&mysolver, factor_threshold);
            time_factor = GetCurrentTime() - tt_factor;
            printf("---- time of factor: %12.6lf ms\n", time_factor);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                tt_solve = GetCurrentTime();
                KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);
                time_solve = GetCurrentTime() - tt_solve;
                printf("---- time of solve: %12.6lf ms\n", time_solve);

#ifdef KML_DSS_IR_
                // IR
                for (int index = 0; index < n; ++index)
                {
                    x[index] = creal(kml_dss_solver_x[index]); // updating solution
                    x_im[index] = cimag(kml_dss_solver_x[index]);
                }
                spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im); // Ax
                for (int index = 0; index < n; ++index)
                {
                    solver_r_re[index] = b[index] - solver_r_re[index]; // solver_r = b - Ax
                    solver_r_im[index] = b_im[index] - solver_r_im[index];
                }
                vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                printf("\n>>>> before IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                while (solver_r_l2_norm / solver_b_l2_norm >= sys_rtol && ir_times < IR_times)
                {
                    ++ir_times;
                    printf("\n>>>> In IR times = %d\n, information", ir_times);

                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_r[index] = solver_r_re[index] + solver_r_im[index] * _Complex_I;
                        kml_dss_solver_e[index] = solver_e_re[index] + solver_e_im[index] * _Complex_I;
                    }
                    KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_r);
                    KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_e);
                    KMLComplexSolverSolve(&mysolver, kml_refine_method, refine_maxit, refine_tol);

                    // updating solution
                    for (int index = 0; index < n; ++index)
                    {
                        kml_dss_solver_x[index] += kml_dss_solver_e[index];
                        x[index] = creal(kml_dss_solver_x[index]);
                        x_im[index] = cimag(kml_dss_solver_x[index]);
                    }

                    spmv_complex(n, row_ptr, col_idx, val, val_im, x, x_im, solver_r_re, solver_r_im);
                    for (int index = 0; index < n; ++index)
                    {
                        solver_r_re[index] = b[index] - solver_r_re[index];
                        solver_r_im[index] = b_im[index] - solver_r_im[index];
                    }
                    vec2norm_complex(solver_r_re, solver_r_im, &solver_r_l2_norm, &solver_r_l2_norm_i, n);
                    printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                    printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                    printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);
                }
                printf("\n>>>> after IR, information:\n");
                printf(">>>> L2 rhs norm: %021.16le\n", solver_b_l2_norm);
                printf(">>>> L2 residual norm: %021.16le\n", solver_r_l2_norm);
                printf(">>>> L2 relative residual norm: %021.16le\n", solver_r_l2_norm / solver_b_l2_norm);

                KMLComplexSolverRHSCreate(&mysolver, n, kml_dss_solver_b);
                KMLComplexSolverSOLCreate(&mysolver, n, kml_dss_solver_x);
                KmlSolverMatrixSetValue(mysolver.dss_solver_b, kml_dss_solver_b);
                KmlSolverMatrixSetValue(mysolver.dss_solver_x, kml_dss_solver_x);
#endif // KML_DSS_IR_
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

        // free ir
        free(solver_r_re);
        free(solver_r_im);
        free(solver_e_re);
        free(solver_e_im);
        free(kml_dss_solver_r);
        free(kml_dss_solver_e);

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

/*
 * command line
 *     -file_mat            <path/to/mat>
 *     -file_rhs            <path/to/rhs>
 *     -type                <0/1/2>
 *     -sys_type            <0/1>
 *     -test_frequency      <int>
 *     -ir_times            <int>
 *     -sys_rtol            <double>
 *     -n_thread            <int>
 *     -nt_rdr              <int>
 *     -factor_threshold    <double>
 *     -matrix_type         <string>
 */
