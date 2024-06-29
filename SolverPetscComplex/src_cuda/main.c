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
    double tt, time;

    int type = 0; // type to output time, 0: end to end time; 1:solver time + solve time; 2:solve time; default 0
    // int test_frequency = 10;  // run code frequency
    int test_frequency = 1; // run code frequency
    int sys_type = 0;       // type of algebraic systems, 0: real, 1: complex; default 0

    int myrank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

#if 0
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-sys_type", argv[index]))
        {
            sys_type = atoi(argv[index + 1]);
        }
        if (strstr("-type", argv[index]))
        {
            type = atoi(argv[index + 1]);
        }
        if (strstr("-test_frequency", argv[index]))
        {
            test_frequency = atoi(argv[index + 1]);
        }
    }
#endif

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    PetscBool sys_flag;
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-sys_type", &sys_type, &sys_flag));
    if (sys_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "value for sys_type: %d\n", sys_type));
    }
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-type", &type, &sys_flag));
    if (sys_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "value for type: %d\n", type));
    }
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-test_frequency", &test_frequency, &sys_flag));
    if (sys_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "value for test_frequency: %d\n", test_frequency));
    }

    char path_mat[PETSC_MAX_PATH_LEN];
    char path_rhs[PETSC_MAX_PATH_LEN];
    PetscBool path_flag;

    PetscCall(PetscOptionsGetString(NULL, NULL, "-file_mat", path_mat, sizeof(path_mat), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Matrix file: %s\n", path_mat));
    }

    PetscCall(PetscOptionsGetString(NULL, NULL, "-file_rhs", path_rhs, sizeof(path_rhs), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "RHS file: %s\n", path_rhs));
    }

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
#endif // DIRECT_SOLVER

    // file process
    Matrix mat_data;
    Vector rhs_data;
    MatrixProcess(path_mat, &mat_data);
    VectorProcess(path_rhs, &rhs_data);

#ifdef ITERATIVE_SOLVER
    if (sys_type == 0) // system
    {
        // MPI iterative solver sample
        MySolver mysolver;
        SolverPetscInitialize(argc, argv, &mat_data, &rhs_data, &mysolver);

        if (type == 0) // check end to end time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                SolverPetscPreprocess(argc, argv, &mysolver);
                SolverPetscSolve(argc, argv, &mysolver);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
            SolverPetscResidualCheck(argc, argv, &mysolver);
        }
        else if (type == 1) // check preprocess + iterative_solver time
        {
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }
        else
        { // check iterative_solver time
            tt = GetCurrentTime();
            for (int i = 0; i < test_frequency; i++)
            {
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            time = (GetCurrentTime() - tt) / (double)(test_frequency);
        }

        PetscCall(KSPDestroy(&(mysolver.ksp)));
        PetscCall(MatDestroy(&(mysolver.solver_a)));
        PetscCall(VecDestroy(&(mysolver.solver_b)));
        PetscCall(VecDestroy(&(mysolver.solver_r)));
        PetscCall(VecDestroy(&(mysolver.solver_x)));
    }
#endif // ITERATIVE_SOLVER

    // free memory
    free(mat_data.row_idx);
    free(mat_data.col_idx);
    free(mat_data.val_re);
    free(mat_data.val_im);
    free(rhs_data.val_re);
    free(rhs_data.val_im);

    /* ========================================== */
    // Step 3: Check time, memory and correctness
    /* ========================================== */
    // check the memory
    // mem_usage();

    if (myrank == 0)
    {

        fprintf(stdout, "------------------------------------------\n");
        if (type == 0)
        {
            fprintf(stdout, "CHECK end to end time :         %12.6lf ms\n", time);
        }
    }

    PetscCall(PetscFinalize());
    MPI_Finalize();
    return 0;
}