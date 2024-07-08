#include "mysolver.h"

void SolverPetscInitialize(int argc, char **argv, MySolver *mysolver)
{
    PetscMPIInt myrank, mysize;
    PetscCallMPI(MPI_Comm_rank(MPI_COMM_WORLD, &myrank));
    PetscCallMPI(MPI_Comm_size(MPI_COMM_WORLD, &mysize));

#if 0
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
#endif

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

    // size of matrix
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "==== basic information of linear system ====\n"));
    int m_mat = 0, n_mat = 0, nnz_mat = 0;
    PetscCall(MatGetSize(mysolver->solver_a, &m_mat, &n_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix Row = %d, matrix Column = %d\n", m_mat, n_mat));
    MatInfo info_mat;
    PetscCall(MatGetInfo(mysolver->solver_a, MAT_GLOBAL_SUM, &info_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "matrix nz_allocated = %ld, matrix nz_used = %ld, matrix nz_unneeded = %ld\n",
                          (long)(info_mat.nz_allocated), (long)(info_mat.nz_used), (long)(info_mat.nz_unneeded)));
    int n_vec = 0;
    PetscCall(VecGetSize(mysolver->solver_b, &n_vec));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "vector Row = %d\n", n_vec));

    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mysolver->ksp)));
}

void SolverPetscPreprocess(int argc, char **argv, MySolver *mysolver)
{
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));
    PetscCall(KSPSetFromOptions(mysolver->ksp));

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
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
}

void SolverPetscResidualCheck(int argc, char **argv, MySolver *mysolver)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "            || b ||_2 = %021.16le\n", b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "|| r ||_2 / || b ||_2 = %021.16le\n", r_norm_2 / b_norm_2));
}

void SolverPetscSolutionFileIO(MySolver *mysolver)
{
    PetscViewer viewer;
    PetscInt n_solver_x = 0;
    const PetscScalar *array_solver_x;
    PetscCall(VecGetSize(mysolver->solver_x, &n_solver_x));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, "answer_x.txt", &viewer));
    PetscCall(PetscViewerASCIIPrintf(viewer, "%d\n", n_solver_x));
    PetscCall(VecGetArrayRead(mysolver->solver_x, &array_solver_x));
    for (int index = 0; index < n_solver_x; ++index)
    {
        PetscCall(PetscViewerASCIIPrintf(viewer, "%021.16le\t%021.16le\n",
                                         PetscRealPart(array_solver_x[index]),
                                         PetscImaginaryPart(array_solver_x[index])));
    }
    PetscCall(VecRestoreArrayRead(mysolver->solver_x, &array_solver_x));
    PetscCall(PetscViewerDestroy(&viewer));
}

void SolverPetscDestroy(MySolver *mysolver)
{
    PetscCall(KSPDestroy(&(mysolver->ksp)));
    PetscCall(MatDestroy(&(mysolver->solver_a)));
    PetscCall(VecDestroy(&(mysolver->solver_b)));
    PetscCall(VecDestroy(&(mysolver->solver_r)));
    PetscCall(VecDestroy(&(mysolver->solver_x)));
}

#if 0
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
#endif
