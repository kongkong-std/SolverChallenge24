#include "mysolver.h"

void SolverPetscInitialize(int argc, char **argv, MySolver *mysolver)
{
    PetscMPIInt irank, nproc;

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &irank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &nproc));

    char *path_mat = NULL, *path_rhs = NULL;

    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-file_mat", argv[index]))
        {
            path_mat = argv[index + 1];
        }
        if (strstr("-file_rhs", argv[index]))
        {
            path_rhs = argv[index + 1];
        }
    }

    PetscViewer fd;

    // matrix file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_mat, FILE_MODE_READ, &fd));
    PetscCall(MatCreate(PETSC_COMM_WORLD, &(mysolver->solver_a)));
    PetscCall(MatLoad(mysolver->solver_a, fd));
    PetscCall(PetscViewerDestroy(&fd));

    // rhs file
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, path_rhs, FILE_MODE_READ, &fd));
    PetscCall(VecCreate(PETSC_COMM_WORLD, &(mysolver->solver_b)));
    PetscCall(VecLoad(mysolver->solver_b, fd));
    PetscCall(PetscViewerDestroy(&fd));

    // sol vector and residual vector
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_x)));
    PetscCall(VecDuplicate(mysolver->solver_b, &(mysolver->solver_r)));
}

void SolverPetscPreprocess(int argc, char **argv, MySolver *mysolver)
{
    PetscMPIInt irank, nproc;

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &irank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &nproc));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mysolver->ksp)));
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));

    // pcshell
#if 0
    PetscCall(KSPGetPC(mysolver->ksp, &(mysolver->pc)));
    PetscBool def_pc_agmg = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-def_pc_agmg", &def_pc_agmg, NULL));
    if (def_pc_agmg)
    {
        PetscCall(PCSetType(mysolver->pc, PCSHELL));

        PetscCall(PCShellSetSetUp(mysolver->pc, AGMGShellPCSetup));

        PetscCall(PCShellSetApply(mysolver->pc, AGMGShellPCApply));
        PetscCall(PCShellSetContext(mysolver->pc, NULL));
    }
#endif
}

void SolverPetscSolve(int argc, char **argv, MySolver *mysolver)
{
    PetscMPIInt irank, nproc;

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &irank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &nproc));

    PetscCall(KSPSetFromOptions(mysolver->ksp));
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
}

void SolverPetscResidualCheck(int argc, char **argv, MySolver *mysolver)
{
    PetscMPIInt irank, nproc;

    PetscFunctionBeginUser;

    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &irank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &nproc));

    PetscReal b_norm_2 = 0., b_norm_1 = 0., b_norm_infty = 0.;
    PetscReal r_norm_2 = 0., r_norm_1 = 0., r_norm_infty = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));
#if 0
    PetscCall(VecNorm(mysolver->solver_b, NORM_1, &b_norm_1));
    PetscCall(VecNorm(mysolver->solver_b, NORM_INFINITY, &b_norm_infty));
#endif

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));
#if 0
    PetscCall(VecNorm(mysolver->solver_r, NORM_1, &r_norm_1));
    PetscCall(VecNorm(mysolver->solver_r, NORM_INFINITY, &r_norm_infty));
#endif

    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "L1-norm: \t|| r || / || b || = %021.16le\n", r_norm_1 / b_norm_1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "L2-norm: \t|| r || / || b || = %021.16le\n \
    \t|| b || _ 2\n",
                          r_norm_2 / b_norm_2, b_norm_2));
    // PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Linfty-norm: \t|| r || / || b || = %021.16le\n", r_norm_infty / b_norm_infty));

#if 0
    PetscViewer fd;                                   /* viewer */
    char file_x[PETSC_MAX_PATH_LEN] = "solution.txt"; /* name of output file with solution vector */
    // ierr = PetscOptionsGetString(NULL, NULL, "-f_x", file_x, sizeof(file_x), NULL); CHKERRQ(ierr);
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_x, &fd));
    PetscCall(PetscViewerPushFormat(fd, PETSC_VIEWER_DEFAULT));
    PetscCall(VecView(mysolver->solver_x, fd));
    PetscCall(PetscViewerPopFormat(fd));
    PetscCall(PetscViewerDestroy(&fd));
#endif

    PetscCall(KSPDestroy(&(mysolver->ksp)));
    PetscCall(MatDestroy(&(mysolver->solver_a)));
    PetscCall(VecDestroy(&(mysolver->solver_b)));
    PetscCall(VecDestroy(&(mysolver->solver_r)));
    PetscCall(VecDestroy(&(mysolver->solver_x)));
}

void SolverPetscGetLinearSystem(const MySolver *mysolver, int *m, int *n, int *nnz,
                                int **row_ptr, int **col_idx, double **val, double **x, double **b)
{
    Mat solver_mat = mysolver->solver_a;
    Vec solver_rhs = mysolver->solver_b, solver_sol = mysolver->solver_x;

    int m_size = 0, n_size = 0, nnz_size = 0;
    int n_loc = 0, nnz_loc = 0;

    PetscCall(MatGetSize(solver_mat, &m_size, &n_size));
    *m = m_size;
    *n = n_size;
    printf("linear system: \t m = %d, n = %d\n", m_size, n_size);

    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    int loc_row_start = 0, loc_row_end = 0;
    PetscBool done = PETSC_TRUE;

    PetscCall(MatGetRowIJ(solver_mat, 0, PETSC_FALSE, PETSC_TRUE, &n_loc,
                          &csr_ia, &csr_ja, &done));
    PetscCall(MatGetOwnershipRange(solver_mat, &loc_row_start, &loc_row_end));

    nnz_loc = csr_ia[n_loc];
    *nnz = nnz_loc;
    printf("\t nnz = %d\n", nnz_loc);
    printf("\t loc_row_start = %d\n", loc_row_start);
    printf("\t loc_row_end = %d\n", loc_row_end);

    int *loc_row_idx = NULL;

    if ((*val = (double *)malloc(nnz_loc * sizeof(double))) == NULL ||
        (*b = (double *)malloc(n_loc * sizeof(double))) == NULL ||
        (*x = (double *)malloc(n_loc * sizeof(double))) == NULL ||
        (*col_idx = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (*row_ptr = (int *)malloc((n_loc + 1) * sizeof(int))) == NULL ||
        (loc_row_idx = (int *)malloc(n_loc * sizeof(int))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'solver linear system get failed\'\n");
        exit(EXIT_FAILURE);
    }

    // assigning value
    /*
     * loc_row_idx = loc_row_restart : loc_row_end - 1
     * row_ptr = csr_ia[0] : csr_ia[n_loc]
     * col_idx = csr_ja[0] : csr_ja[nnz_loc - 1]
     * val = csr format assigning value
     * b = solver_rhs[loc_row_start] : solver_rhs[loc_row_end - 1]
     * x = solver_sol[loc_row_start] : solver_sol[loc_row_end - 1]
     */
    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        loc_row_idx[index - loc_row_start] = index;
    }
    for (int index = 0; index < n_loc + 1; ++index)
    {
        (*row_ptr)[index] = csr_ia[index];
    }
    for (int index = 0; index < nnz_loc; ++index)
    {
        (*col_idx)[index] = csr_ja[index];
    }
    for (int index = 0; index < n_loc; ++index)
    {
        // matrix element
        int index_start = csr_ia[index];
        int index_end = csr_ia[index + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(solver_mat, loc_row_idx[index], csr_ja[index_j], &val_tmp));
            (*val)[index_j] = val_tmp;
        }

        // rhs element
        PetscScalar val_tmp1;
        PetscCall(VecGetValues(solver_rhs, 1, loc_row_idx + index, &val_tmp1));
        (*b)[index] = val_tmp1;

        // solution element
        PetscScalar val_tmp2;
        PetscCall(VecGetValues(solver_sol, 1, loc_row_idx + index, &val_tmp2));
        (*x)[index] = val_tmp2;
    }
}

//! real system
void analyse(MySolver *solver, const int n, const int *row_ptr, const int *col_idx)
{
}

void preprocess(MySolver *solver, const int n, const double *val)
{
}

void iterative_solver(MySolver *solver, const int n, const double *x, const double *b)
{
}

//! complex system
void analyse_complex(MySolverComplex *solver, const int n, const int *row_ptr, const int *col_idx)
{
}

void preprocess_complex(MySolverComplex *solver, const int n, const double *val, const double *val_im)
{
}

void iterative_solver_complex(MySolverComplex *solver, const int n, const double *x,
                              const double *x_im, const double *b, const double *b_im)
{
}
