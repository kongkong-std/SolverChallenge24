// system file
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/time.h>

// read matrix file
#include "mmio_highlevel.h"

#include <mpi.h>

// utlise function file
#include "utlise.h"
#include "utlise_long.h"

#include "mysolver.h"

int main(int argc, char **argv)
{
    int m, n, nnzA, isSymmetricA;
    int *row_ptr = NULL; // the csr row pointer array of matrix A
    int *col_idx = NULL; // the csr column index array of matrix A

    double *val = NULL;    // the csr value array of matrix A (real number)
    double *val_im = NULL; // the csr value array of matrix A (imaginary number)

    double *x = NULL;
    double *x_im = NULL; // solution vector x, (x: real number, x_im: imaginary number)

    double *b = NULL;
    double *b_im = NULL; // right-hand side vector b, (b: real number, b_im: imaginary number)

    double tt, time;

    char *filename_matrix;    // the filename of matrix A
    char *filename_b;         // the filename of right-hand side vector b
    int read_matrix_base = 1; // 0-base or 1-base, default 1-base
    int type = 0;             // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    int test_frequency = 10;  // run code frequency
    int sys_type = 0;         // type of algebraic systems, 0: real, 1: complex; default 0

    int myrank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for(int index = 0; index < argc; ++index)
    {
        if(strstr("-sys_type", argv[index]))
        {
            sys_type = atoi(argv[index + 1]);
        }
    }

#if 0
    /* ========================================== */
    // Step 0: Read command line argument
    /* ========================================== */
    if (argc == 3)
    {
        filename_matrix = argv[1];
        filename_b = argv[2];
    }
    else if (argc == 4)
    {
        filename_matrix = argv[1];
        filename_b = argv[2];
        read_matrix_base = atoi(argv[3]);

        if (!((read_matrix_base == 0) || (read_matrix_base == 1)))
        {
            if (myrank == 0)
            {
                fprintf(stdout, "read_base error input %d\n", read_matrix_base);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
    }
    else if (argc == 5)
    {
        filename_matrix = argv[1];
        filename_b = argv[2];
        read_matrix_base = atoi(argv[3]);
        if (!((read_matrix_base == 0) || (read_matrix_base == 1)))
        {
            if (myrank == 0)
            {
                fprintf(stdout, "read_base error input %d\n", read_matrix_base);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
        type = atoi(argv[4]);
        if (type > 2 || type < 0)
        {
            if (myrank == 0)
            {
                fprintf(stdout, "type error input %d\n", type);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
    }
    else if (argc == 6)
    {
        filename_matrix = argv[1];
        filename_b = argv[2];
        read_matrix_base = atoi(argv[3]);
        if (!((read_matrix_base == 0) || (read_matrix_base == 1)))
        {
            if (myrank == 0)
            {
                fprintf(stdout, "read_base error input %d\n", read_matrix_base);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
        type = atoi(argv[4]);
        if (type > 2 || type < 0)
        {
            if (myrank == 0)
            {
                fprintf(stdout, "type error input %d\n", type);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
        sys_type = atoi(argv[5]);
        if (sys_type > 1 || sys_type < 0)
        {
            if (myrank == 0)
            {
                fprintf(stdout, "sys_type error input %d\n", sys_type);
                print_help();
            }
            MPI_Finalize();
            return -1;
        }
    }
    else
    {
        if (myrank == 0)
        {
            print_help();
        }
        MPI_Finalize();
        return -1;
    }

    //! If the user does not specify parameter sys_type,
    //! it will automatically determine whether it is a real system or a complex system according to the mtx file
    if (argc < 6)
        read_mtx_header(filename_matrix, &sys_type);

    if (myrank == 0)
        fprintf(stdout, "matrix name :      %s\nvectorb name :     %s\nread function :    base-%d\ntype :             %d\nsys_type :         %d\n",
                filename_matrix, filename_b, read_matrix_base, type, sys_type);
#endif

#if 0
    /* ========================================== */
    // Step 1: Load matrix and rhs from mtx files
    /* ========================================== */
    if (myrank == 0)
    {
        if (sys_type == 0) // real system
        {
            // load matrix
            mmio_allinone(&m, &n, &nnzA, &isSymmetricA, &read_matrix_base, &row_ptr, &col_idx, &val, filename_matrix);
            if (m != n)
            {
                fprintf(stdout, "Invalid matrix size.\n");
                MPI_Finalize();
                return 0;
            }

            x = (double *)malloc(sizeof(double) * n);
            b = (double *)malloc(sizeof(double) * n);

            // load right-hand side vector b
            load_vector(n, b, filename_b);
        }
        else
        { // complex system
            mmio_allinone_complex(&m, &n, &nnzA, &isSymmetricA, &read_matrix_base, &row_ptr, &col_idx, &val, &val_im, filename_matrix);
            if (m != n)
            {
                fprintf(stdout, "Invalid matrix size.\n");
                return 0;
            }

            x = (double *)malloc(sizeof(double) * n);
            x_im = (double *)malloc(sizeof(double) * n);
            b = (double *)malloc(sizeof(double) * n);
            b_im = (double *)malloc(sizeof(double) * n);

            // load right-hand side vector b
            load_b_complex(n, b, b_im, filename_b);
        }
    }
#endif

/* ========================================== */
// Step 2: Solve the linear system
/* ========================================== */
#ifdef DIRECT_SOLVER
    if (sys_type == 0) // real system
    {
        // MPI direct solver sample
        MySolver mysolver;

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                preprocess(&mysolver, n, row_ptr, col_idx);
                direct_solver(&mysolver, n, val);
                solve(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            preprocess(&mysolver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                direct_solver(&mysolver, n, val);
                solve(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time

            preprocess(&mysolver, n, row_ptr, col_idx);
            direct_solver(&mysolver, n, val);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                solve(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
    }
    else
    { // complex system
      // MPI direct solver sample
        MySolverComplex mysolver;

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                preprocess_complex(&mysolver, n, row_ptr, col_idx);
                direct_solver_complex(&mysolver, n, val, val_im);
                solve_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check direct_solver + solve time
        {
            preprocess_complex(&mysolver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                direct_solver_complex(&mysolver, n, val, val_im);
                solve_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check solve time

            preprocess_complex(&mysolver, n, row_ptr, col_idx);
            direct_solver_complex(&mysolver, n, val, val_im);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);

                solve_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
    }
#endif

#ifdef ITERATIVE_SOLVER
    if (sys_type == 0) // real system
    {
        // MPI iterative solver sample
        MySolver mysolver;

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                }
                MPI_Barrier(MPI_COMM_WORLD);

                analyse(&mysolver, n, row_ptr, col_idx);
                preprocess(&mysolver, n, val);
                iterative_solver(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check preprocess + iterative_solver time
        {
            analyse(&mysolver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                }
                MPI_Barrier(MPI_COMM_WORLD);

                preprocess(&mysolver, n, val);
                iterative_solver(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check iterative_solver time

            analyse(&mysolver, n, row_ptr, col_idx);
            preprocess(&mysolver, n, val);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                }
                MPI_Barrier(MPI_COMM_WORLD);

                iterative_solver(&mysolver, n, x, b);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
    }
    else
    { // complex system
      // MPI iterative solver sample
        MySolverComplex mysolver;

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                    memset(x_im, 0.0, sizeof(double) * n);
                }
                MPI_Barrier(MPI_COMM_WORLD);

                analyse_complex(&mysolver, n, row_ptr, col_idx);
                preprocess_complex(&mysolver, n, val, val_im);
                iterative_solver_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else if (type == 1) // check preprocess + iterative_solver time
        {
            analyse_complex(&mysolver, n, row_ptr, col_idx);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                    memset(x_im, 0.0, sizeof(double) * n);
                }
                MPI_Barrier(MPI_COMM_WORLD);

                preprocess_complex(&mysolver, n, val, val_im);
                iterative_solver_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check iterative_solver time

            analyse_complex(&mysolver, n, row_ptr, col_idx);
            preprocess_complex(&mysolver, n, val, val_im);

            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                if (myrank == 0)
                {
                    // Note that calling the iterative method in a loop requires initializing the solution vector to 0 each time
                    memset(x, 0.0, sizeof(double) * n); // initial vector x
                    memset(x_im, 0.0, sizeof(double) * n);
                }
                MPI_Barrier(MPI_COMM_WORLD);

                iterative_solver_complex(&mysolver, n, x, x_im, b, b_im);

                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
    }
#endif

    /* ========================================== */
    // Step 3: Check time, memory and correctness
    /* ========================================== */
    // check the memory
    mem_usage();

    if (myrank == 0)
    {

        fprintf(stdout, "------------------------------------------\n");

        // check the memory
        // mem_usage();

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

        // store x to a file
        char *answer_x = "answer_x.rhs";
        if (sys_type == 0) // real rhs
            store_x(n, x, answer_x);
        else // complex rhs
            store_x_complex(n, x, x_im, answer_x);

        // check the correctness
        if (sys_type == 0) // real system
        {
            // 1) using double precision
            check_correctness(n, row_ptr, col_idx, val, x, b);
            // 2) using long double precision
            check_correctness_ld_d2ld(n, row_ptr, col_idx, val, x, b);
        }
        else
        { // complex system
            // 1) using double precision
            check_correctness_complex(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);
            // 2) using long double precision
            check_correctness_complex_ld_d2ld(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);

            free(val_im);
            free(x_im);
            free(b_im);
        }
    }

    if (myrank == 0)
    {
        free(row_ptr);
        free(col_idx);
        free(val);
        free(x);
        free(b);
    }

    MPI_Finalize();
    return 0;
}