// SPDX-FileCopyrightText: 2011 - 2024 NVIDIA CORPORATION. All Rights Reserved.
//
// SPDX-License-Identifier: BSD-3-Clause

#define MAX_MSG_LEN 4096

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// #include "cuda_runtime.h"

/* CUDA error macro */
#define CUDA_SAFE_CALL(call)                                              \
    do                                                                    \
    {                                                                     \
        cudaError_t err = call;                                           \
        if (cudaSuccess != err)                                           \
        {                                                                 \
            fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err));         \
            exit(EXIT_FAILURE);                                           \
        }                                                                 \
    } while (0)

/* standard or dynamically load library */
#ifdef AMGX_DYNAMIC_LOADING
#include "amgx_capi.h"
#else
#include "amgx_c.h"
#endif

#include "utlise.h"
#include "file_process.h"

const char *GetAMGXErrorString(AMGX_RC rc);

/* print error message and exit */
void errAndExit(const char *err)
{
    printf("%s\n", err);
    fflush(stdout);
    exit(1);
}

/* print callback (could be customized) */
void print_callback(const char *msg, int length)
{
    printf("%s", msg);
}

/* print usage and exit */
void printUsageAndExit()
{
    printf("%s", "Usage: ./amgx_capi [-mode [hDDI | hDFI | hFFI | dDDI | dDFI | dFFI]] [-m file] [-c config_file] [-amg \"variable1=value1 ... variable3=value3\"]\n");
    printf("%s", "     -mode:   select the solver mode\n");
    printf("%s", "     -m file: read matrix stored in the file\n");
    printf("%s", "     -c:      set the amg solver options from the config file\n");
    printf("%s", "     -amg:    set the amg solver options from the command line\n");
    exit(0);
}

/* parse parameters */
int findParamIndex(const char **argv, int argc, const char *parm)
{
    int count = 0;
    int index = -1;

    for (int i = 0; i < argc; i++)
    {
        if (strncmp(argv[i], parm, 100) == 0)
        {
            index = i;
            count++;
        }
    }

    if (count == 0 || count == 1)
    {
        return index;
    }
    else
    {
        printf("Error, parameter %s has been specified more than once, exiting\n", parm);
        exit(1);
    }

    return -1;
}

/* reade geometry (advanced input parmeters) */
void readGeometry(const char *fname, double **geo_x, double **geo_y, double **geo_z, int *dim, int *numrows)
{
    printf("Reading geometry from file: '%s'\n", fname);
    FILE *fin = fopen(fname, "r");

    if (!fin)
    {
        printf("Error opening file '%s'\n", fname);
        exit(1);
    }

    int n, dimension;

    if (2 != fscanf(fin, "%d %d\n", &n, &dimension))
    {
        errAndExit("Bad format\n");
    }

    *geo_x = (double *)malloc(n * sizeof(double));
    *geo_y = (double *)malloc(n * sizeof(double));

    if (dimension == 3)
    {
        *geo_y = (double *)malloc(n * sizeof(double));

        for (int i = 0; i < n; i++)
            if (3 != fscanf(fin, "%lf %lf %lf\n", *geo_x + i, *geo_y + i, *geo_z + i))
            {
                errAndExit("Bad format\n");
            }
    }
    else if (dimension == 2)
    {
        for (int i = 0; i < n; i++)
            if (2 != fscanf(fin, "%lf %lf\n", *geo_x + i, *geo_y + i))
            {
                errAndExit("Bad format\n");
            }
    }

    *dim = dimension;
    *numrows = n;
}

/* reade coloring (advanced input parmeters) */
void readColoring(const char *fname, int **row_coloring, int *colored_rows, int *num_colors)
{
    printf("Reading coloring from file: '%s'\n", fname);
    FILE *fin = fopen(fname, "r");

    if (!fin)
    {
        printf("Error opening file '%s'\n", fname);
        exit(1);
    }

    int n, colors_num;

    if (2 != fscanf(fin, "%d %d\n", &n, &colors_num))
    {
        errAndExit("Bad format\n");
    }

    *row_coloring = (int *)malloc(n * sizeof(int));

    for (int i = 0; i < n; i++)
        if (1 != fscanf(fin, "%d\n", *row_coloring + i))
        {
            errAndExit("Bad format\n");
        }

    *colored_rows = n;
    *num_colors = colors_num;
}

int main(int argc, const char **argv)
{
    // parameter parsing
    int pidx = 0;
    int pidy = 0;

    // versions
    int major, minor;
    char *ver, *date, *time;

    // file path to matrix and vector
    char *path_mat = NULL;
    char *path_rhs_t0 = NULL, *path_rhs_t1 = NULL, *path_rhs_t2 = NULL;

    // input geometry
    double *gx = NULL;
    double *gy = NULL;
    double *gz = NULL;

    // input coloring
    int dim = 0;
    int numrows = 0;
    int num_colors = 0;
    int colored_rows = 0;
    int *row_coloring = NULL;

    // input matrix and rhs/solution
    int bsize_x = 0;
    int bsize_y = 0;
    int sol_size = 0;
    int sol_bsize = 0;
    int m = 0, n = 0, nnz = 0; // size and number of non-zeros
#if 1
    double tt_amgx = 0;
    int *row_ptr = NULL, *col_idx = NULL; // row offsets and column indices
    double *val = NULL;                   // non-zero values

    double *rhs_t0 = NULL, *rhs_t1 = NULL, *rhs_t2 = NULL; // rhs of t0, t1, t2

    // file process
    if ((pidx = findParamIndex(argv, argc, "-path_matrix")) != -1)
    {
        path_mat = argv[pidx + 1];
        RealCOO2CSRMatrixFileProcess(path_mat, &m, &n, &nnz, &row_ptr, &col_idx, &val);
        if (m != n)
        {
            fprintf(stderr, "Invalid matrix size.\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        errAndExit("ERROR: no matrix was specified");
    }

    if ((pidx = findParamIndex(argv, argc, "-path_rhs_t0")) != -1)
    {
        path_rhs_t0 = argv[pidx + 1];
        RealRHSFileProcess(path_rhs_t0, &rhs_t0);
    }
    else
    {
        errAndExit("ERROR: no rhs_t0 was specified");
    }

    if ((pidx = findParamIndex(argv, argc, "-path_rhs_t1")) != -1)
    {
        path_rhs_t1 = argv[pidx + 1];
        RealRHSFileProcess(path_rhs_t1, &rhs_t1);
    }
    else
    {
        errAndExit("ERROR: no rhs_t1 was specified");
    }

    if ((pidx = findParamIndex(argv, argc, "-path_rhs_t2")) != -1)
    {
        path_rhs_t2 = argv[pidx + 1];
        RealRHSFileProcess(path_rhs_t2, &rhs_t2);
    }
    else
    {
        errAndExit("ERROR: no rhs_t2 was specified");
    }
#endif // matrix file and 3 rhs file process

    // library handles
    AMGX_Mode mode;
    AMGX_config_handle cfg;
    AMGX_resources_handle rsrc;
    AMGX_matrix_handle A;
    AMGX_vector_handle b, x;
#if 1
    AMGX_vector_handle b_0, x_0; // t0
    AMGX_vector_handle b_1, x_1; // t1
    AMGX_vector_handle b_2, x_2; // t2
#endif                           // 3 rhs and 3 solution
    AMGX_solver_handle solver;

    // status handling
    AMGX_SOLVE_STATUS status;

    /* check arguments */
    if (argc == 1)
    {
        printUsageAndExit();
    }

    /* load the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    void *lib_handle = NULL;
    // open the library
#ifdef _WIN32
    lib_handle = amgx_libopen("amgxsh.dll");
#else
    lib_handle = amgx_libopen("libamgxsh.so");
#endif

    if (lib_handle == NULL)
    {
        errAndExit("ERROR: can not load the library");
    }

    // load all the routines
    if (amgx_liblink_all(lib_handle) == 0)
    {
        amgx_libclose(lib_handle);
        errAndExit("ERROR: corrupted library loaded\n");
    }

#endif
    /* init */
    AMGX_SAFE_CALL(AMGX_initialize());
    /* system */
    AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
    AMGX_SAFE_CALL(AMGX_install_signal_handler());

    /* get api and build info */
    if ((pidx = findParamIndex(argv, argc, "--version")) != -1)
    {
        AMGX_get_api_version(&major, &minor);
        printf("amgx api version: %d.%d\n", major, minor);
        AMGX_get_build_info_strings(&ver, &date, &time);
        printf("amgx build version: %s\nBuild date and time: %s %s\n", ver, date, time);
        AMGX_SAFE_CALL(AMGX_finalize());
        /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
        amgx_libclose(lib_handle);
#endif
        exit(0);
    }

    /* get mode */
    if ((pidx = findParamIndex(argv, argc, "-mode")) != -1)
    {
        if (strncmp(argv[pidx + 1], "hDDI", 100) == 0)
        {
            mode = AMGX_mode_hDDI;
        }
        else if (strncmp(argv[pidx + 1], "hDFI", 100) == 0)
        {
            mode = AMGX_mode_hDFI;
        }
        else if (strncmp(argv[pidx + 1], "hFFI", 100) == 0)
        {
            mode = AMGX_mode_hFFI;
        }
        else if (strncmp(argv[pidx + 1], "dDDI", 100) == 0)
        {
            mode = AMGX_mode_dDDI;
        }
        else if (strncmp(argv[pidx + 1], "dDFI", 100) == 0)
        {
            mode = AMGX_mode_dDFI;
        }
        else if (strncmp(argv[pidx + 1], "dFFI", 100) == 0)
        {
            mode = AMGX_mode_dFFI;
        }
        else if (strncmp(argv[pidx + 1], "hCCI", 100) == 0)
        {
            mode = AMGX_mode_hZZI;
        }
        else if (strncmp(argv[pidx + 1], "hZCI", 100) == 0)
        {
            mode = AMGX_mode_hZCI;
        }
        else if (strncmp(argv[pidx + 1], "hZZI", 100) == 0)
        {
            mode = AMGX_mode_hZZI;
        }
        else if (strncmp(argv[pidx + 1], "dCCI", 100) == 0)
        {
            mode = AMGX_mode_dCCI;
        }
        else if (strncmp(argv[pidx + 1], "dZCI", 100) == 0)
        {
            mode = AMGX_mode_dZCI;
        }
        else if (strncmp(argv[pidx + 1], "dZZI", 100) == 0)
        {
            mode = AMGX_mode_dZZI;
        }
        else
        {
            errAndExit("ERROR: invalid mode");
        }
    }
    else
    {
        printf("Warning: No mode specified, using dDDI by default.\n");
        mode = AMGX_mode_dDDI;
    }

    /* create config */
    pidx = findParamIndex(argv, argc, "-amg");
    pidy = findParamIndex(argv, argc, "-c");

    if ((pidx != -1) && (pidy != -1))
    {
        printf("%s\n", argv[pidx + 1]);
        AMGX_SAFE_CALL(AMGX_config_create_from_file_and_string(&cfg, argv[pidy + 1], argv[pidx + 1]));
    }
    else if (pidy != -1)
    {
        AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, argv[pidy + 1]));
    }
    else if (pidx != -1)
    {
        printf("%s\n", argv[pidx + 1]);
        AMGX_SAFE_CALL(AMGX_config_create(&cfg, argv[pidx + 1]));
    }
    else
    {
        errAndExit("ERROR: no config was specified");
    }

    /* example of how to handle errors */
    // char msg[MAX_MSG_LEN];
    // AMGX_RC err_code = AMGX_resources_create_simple(NULL, cfg);
    // AMGX_SAFE_CALL(AMGX_get_error_string(err_code, msg, MAX_MSG_LEN));
    // printf("ERROR: %s\n",msg);
    /* switch on internal error handling (no need to use AMGX_SAFE_CALL after this point) */
    // AMGX_SAFE_CALL(AMGX_config_add_parameters(&cfg, "exception_handling=1"));
    /* create resources, matrix, vector and solver */
    AMGX_resources_create_simple(&rsrc, cfg);
    AMGX_matrix_create(&A, rsrc, mode);
    AMGX_vector_create(&x, rsrc, mode);
    AMGX_vector_create(&b, rsrc, mode);
    AMGX_solver_create(&solver, rsrc, mode, cfg);
#if 1
    AMGX_vector_create(&x_0, rsrc, mode);
    AMGX_vector_create(&b_0, rsrc, mode);
    AMGX_vector_create(&x_1, rsrc, mode);
    AMGX_vector_create(&b_1, rsrc, mode);
    AMGX_vector_create(&x_2, rsrc, mode);
    AMGX_vector_create(&b_2, rsrc, mode);

    AMGX_RC ierr;
    tt_amgx = GetCurrentTime();
    ierr = AMGX_matrix_upload_all(A, n, nnz, 1, 1, row_ptr, col_idx, val, NULL);
    printf(">>>> matrix data >>>> host to device: %12.6lf ms\n", GetCurrentTime() - tt_amgx);
    if (ierr != AMGX_RC_OK)
    {
        fprintf(stderr, "AMGX matrix assigning value failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
        exit(EXIT_FAILURE);
    }

    tt_amgx = GetCurrentTime();
    ierr = AMGX_vector_upload(b_0, n, 1, rhs_t0);
    printf(">>>> rhs t0 data >>>> host to device: %12.6lf ms\n", GetCurrentTime() - tt_amgx);
    if (ierr != AMGX_RC_OK)
    {
        fprintf(stderr, "AMGX rhs_t0 assigning value failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
        exit(EXIT_FAILURE);
    }

    tt_amgx = GetCurrentTime();
    ierr = AMGX_vector_upload(b_1, n, 1, rhs_t1);
    printf(">>>> rhs t1 data >>>> host to device: %12.6lf ms\n", GetCurrentTime() - tt_amgx);
    if (ierr != AMGX_RC_OK)
    {
        fprintf(stderr, "AMGX rhs_t1 assigning value failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
        exit(EXIT_FAILURE);
    }

    tt_amgx = GetCurrentTime();
    ierr = AMGX_vector_upload(b_2, n, 1, rhs_t2);
    printf(">>>> rhs t2 data >>>> host to device: %12.6lf ms\n", GetCurrentTime() - tt_amgx);
    if (ierr != AMGX_RC_OK)
    {
        fprintf(stderr, "AMGX rhs_t2 assigning value failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
        exit(EXIT_FAILURE);
    }

    // initialize x_0 as 0
    AMGX_vector_set_zero(x_0, n, 1);
#endif // 3 rhs and solution create

#if 1
    // free memory
    free(row_ptr);
    free(col_idx);
    free(val);
    free(rhs_t0);
    free(rhs_t1);
    free(rhs_t2);
#endif // 3 rhs and solution

    /* read the input system: matrix [and rhs & solution]
       Please refer to AMGX_read_system description in the AMGX_Reference.pdf
       manual for details on how to specify the rhs and the solution inside
       the input file. If these are not specified than rhs=[1,...,1]^T and
       (initial guess) sol=[0,...,0]^T. */
    if ((pidx = findParamIndex(argv, argc, "-m")) != -1)
    {
        AMGX_read_system(A, b, x, argv[pidx + 1]);
        AMGX_matrix_get_size(A, &n, &bsize_x, &bsize_y);
        AMGX_vector_get_size(x, &sol_size, &sol_bsize);

        if (sol_size == 0 || sol_bsize == 0)
        {
            AMGX_vector_set_zero(x, n, bsize_x);
        }
    }
#if 0
    else
    {
        errAndExit("ERROR: no linear system was specified");
    }
#endif // matrix and 3 rhs file process

    // example of getting initial residual norm
    {
#if 0
        void *t_norm = calloc(bsize_x, sizeof(double));
        AMGX_solver_calculate_residual_norm(solver, A, b, x, t_norm);
        printf("Initial norm: ");
        for (int i = 0; i < bsize_x; i++)
        {
            printf("%12.6e ", ((double *)t_norm)[i]);
        }
        printf("\n");
        free(t_norm);
#endif // matrix and 3 rhs file process

        printf("Initial || rhs ||_2:\n");
        double t_norm = 0.;
        AMGX_solver_calculate_residual_norm(solver, A, b_0, x_0, &t_norm);
        printf(">>>> || b_0 || _ 2 = %021.16le\n", t_norm);

        AMGX_solver_calculate_residual_norm(solver, A, b_1, x_0, &t_norm);
        printf(">>>> || b_1 || _ 2 = %021.16le\n", t_norm);

        AMGX_solver_calculate_residual_norm(solver, A, b_2, x_0, &t_norm);
        printf(">>>> || b_2 || _ 2 = %021.16le\n", t_norm);
    }

    /* read the input geometry */
    if ((pidx = findParamIndex(argv, argc, "-geo")) != -1)
    {
        readGeometry(argv[pidx + 1], &gx, &gy, &gz, &dim, &numrows);

        if (dim != 3)
        {
            gz = NULL;
        }

        if (dim != 0)
        {
            AMGX_matrix_attach_geometry(A, gx, gy, gz, numrows);
        }

        if (gx)
        {
            free(gx);
        }

        if (gy)
        {
            free(gy);
        }

        if (gz)
        {
            free(gz);
        }
    }

    /* read the input coloring */
    if ((pidx = findParamIndex(argv, argc, "-color")) != -1)
    {
        readColoring(argv[pidx + 1], &row_coloring, &colored_rows, &num_colors);

        if (num_colors != 0)
        {
            AMGX_matrix_attach_coloring(A, row_coloring, colored_rows, num_colors);
        }

        if (row_coloring)
        {
            free(row_coloring);
        }
    }

    int rep_times = 2;
    /* read the input coloring */
    if ((pidx = findParamIndex(argv, argc, "-rep")) != -1)
    {
        rep_times = atoi(argv[pidx + 1]);
        printf("rep_times =%d\n", rep_times);
    }

    double tt, elapsed;
    tt = GetCurrentTime();
    for (int sol_time = 0; sol_time < rep_times; ++sol_time)
    {
#if 0
        AMGX_vector_set_zero(x, n, bsize_x);
        if (sol_time > 0)
        {
            AMGX_solver_destroy(solver);
            AMGX_solver_create(&solver, rsrc, mode, cfg);
        }
#endif // matrix and 3 rhs file process

        /* solver setup */
        ierr = AMGX_solver_setup(solver, A);
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX solver setup failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }

        /* solver solve */
        ierr = AMGX_solver_solve(solver, b_0, x_0); // linear system t0
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX solver solve t0 linear system failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }

        double *x_tmp = NULL;
        if ((x_tmp = malloc(n * sizeof(double))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - \'solution tmp vector\'\n");
            exit(EXIT_FAILURE);
        }

        tt_amgx = GetCurrentTime();
        ierr = AMGX_vector_download(x_0, x_tmp);
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX tmp solution data device to host failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }
        printf(">>>> AMGX tmp solution data >>>> device to host %12.6lf ms\n", GetCurrentTime() - tt_amgx);

        tt_amgx = GetCurrentTime();
        ierr = AMGX_vector_upload(x_1, n, 1, x_tmp);
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX tmp solution data host to device failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }
        printf(">>>> AMGX tmp solution data >>>> host to device %12.6lf ms\n", GetCurrentTime() - tt_amgx);

        ierr = AMGX_solver_solve(solver, b_1, x_1); // linear system t1
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX solver solve t1 linear system failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }

        tt_amgx = GetCurrentTime();
        ierr = AMGX_vector_download(x_1, x_tmp);
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX tmp solution data device to host failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }
        printf(">>>> AMGX tmp solution data >>>> device to host %12.6lf ms\n", GetCurrentTime() - tt_amgx);

        tt_amgx = GetCurrentTime();
        ierr = AMGX_vector_upload(x_2, n, 1, x_tmp);
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX tmp solution data host to device failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }
        printf(">>>> AMGX tmp solution data >>>> host to device %12.6lf ms\n", GetCurrentTime() - tt_amgx);

        ierr = AMGX_solver_solve(solver, b_2, x_2); // linear system t2
        if (ierr != AMGX_RC_OK)
        {
            fprintf(stderr, "AMGX solver solve t1 linear system failed! - ERROR: %s\n", GetAMGXErrorString(ierr));
            exit(EXIT_FAILURE);
        }

        free(x_tmp);
    }
    elapsed = (GetCurrentTime() - tt) / (double)(rep_times);
    printf("CHECK end to end time :         %12.6lf ms\n", elapsed);
    /* example of how to change parameters between non-linear iterations */
    // AMGX_config_add_parameters(&cfg, "config_version=2, default:tolerance=1e-12");
    // AMGX_solver_solve(solver, b, x);
    AMGX_solver_get_status(solver, &status);

    printf("n =%d\n", n);
#if 0
    printf("bsize_x =%d\n", bsize_x);
    void *h_result = malloc(n * bsize_x * sizeof(double));
    AMGX_vector_download(x, h_result);
    char *answer_x = "answer_x.rhs";
    store_x(n, h_result, answer_x);
    free(h_result);
#endif // 3 rhs and solution file process
    void *h_result = malloc(n * sizeof(double));
    AMGX_vector_download(x_0, h_result);
    char *answer_x0 = "answer_x_0.rhs";
    store_x(n, h_result, answer_x0);

    AMGX_vector_download(x_1, h_result);
    char *answer_x1 = "answer_x_1.rhs";
    store_x(n, h_result, answer_x1);

    AMGX_vector_download(x_2, h_result);
    char *answer_x2 = "answer_x_02.rhs";
    store_x(n, h_result, answer_x2);
    free(h_result);

    // example of getting initial residual norm
    {
#if 0
        void *t_norm = calloc(bsize_x, sizeof(double));
        AMGX_solver_calculate_residual_norm(solver, A, b, x, t_norm);
        printf("final norm: ");
        for (int i = 0; i < bsize_x; i++)
        {
            printf("%12.6e ", ((double *)t_norm)[i]);
        }
        printf("\n");
        free(t_norm);
#endif // 3 rhs and solution file process
        double t_norm = 0.;
        AMGX_solver_calculate_residual_norm(solver, A, b_0, x_0, &t_norm);
        printf("Final norm:\n");
        printf(">>>> t0 linear system: %021.16le\n", t_norm);

        AMGX_solver_calculate_residual_norm(solver, A, b_1, x_1, &t_norm);
        printf(">>>> t1 linear system: %021.16le\n", t_norm);

        AMGX_solver_calculate_residual_norm(solver, A, b_2, x_2, &t_norm);
        printf(">>>> t2 linear system: %021.16le\n", t_norm);
    }

    /* example of how to print the residual history */
    // int nit;
    // double res;
    // AMGX_solver_get_iterations_number(solver, &nit);
    // for (int i=0; i<nit; i++) {
    //   printf("residual from iteration %d=", i);
    //   for (int j=0; j<bsize_y; j++) {
    //     AMGX_solver_get_iteration_residual(solver, i, j, &res);
    //     printf("%f ", (float)(res));
    //   }
    //   printf("\n");
    // }
    /* example of how to write the linear system to the output */
    // AMGX_write_system(A, b, x, "output.system.mtx");
    /* destroy resources, matrix, vector and solver */
    AMGX_solver_destroy(solver);
    AMGX_vector_destroy(x);
    AMGX_vector_destroy(b);
    AMGX_matrix_destroy(A);
    AMGX_resources_destroy(rsrc);
    /* destroy config (need to use AMGX_SAFE_CALL after this point) */
    AMGX_SAFE_CALL(AMGX_config_destroy(cfg));
    /* shutdown and exit */
    AMGX_SAFE_CALL(AMGX_finalize());
    /* close the library (if it was dynamically loaded) */
#ifdef AMGX_DYNAMIC_LOADING
    amgx_libclose(lib_handle);
#endif
    mem_usage();
    return status;
}

const char *GetAMGXErrorString(AMGX_RC rc)
{
    switch (rc)
    {
    case AMGX_RC_OK:
        return "AMGX_RC_OK";
    case AMGX_RC_BAD_PARAMETERS:
        return "AMGX_RC_BAD_PARAMETERS";
    case AMGX_RC_UNKNOWN:
        return "AMGX_RC_UNKNOWN";
    case AMGX_RC_NOT_SUPPORTED_TARGET:
        return "AMGX_RC_NOT_SUPPORTED_TARGET";
    case AMGX_RC_NOT_SUPPORTED_BLOCKSIZE:
        return "AMGX_RC_NOT_SUPPORTED_BLOCKSIZE";
    case AMGX_RC_CUDA_FAILURE:
        return "AMGX_RC_CUDA_FAILURE";
    case AMGX_RC_THRUST_FAILURE:
        return "AMGX_RC_THRUST_FAILURE";
    case AMGX_RC_NO_MEMORY:
        return "AMGX_RC_NO_MEMORY";
    case AMGX_RC_IO_ERROR:
        return "AMGX_RC_IO_ERROR";
    case AMGX_RC_BAD_MODE:
        return "AMGX_RC_BAD_MODE";
    case AMGX_RC_CORE:
        return "AMGX_RC_CORE";
    case AMGX_RC_PLUGIN:
        return "AMGX_RC_PLUGIN";
    case AMGX_RC_BAD_CONFIGURATION:
        return "AMGX_RC_BAD_CONFIGURATION";
    case AMGX_RC_NOT_IMPLEMENTED:
        return "AMGX_RC_NOT_IMPLEMENTED";
    case AMGX_RC_LICENSE_NOT_FOUND:
        return "AMGX_RC_LICENSE_NOT_FOUND";
    case AMGX_RC_INTERNAL:
        return "AMGX_RC_INTERNAL";
    default:
        return "Unknown error code";
    }
}
