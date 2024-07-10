#include "mysolver.h"

void SolverPetscInitialize(int argc, char **argv, MySolver *mysolver)
{
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
    int myrank, mysize;
    MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    PetscInt n_solver_x = 0;
    PetscInt x_row_start = 0, x_row_end = 0;
    const PetscScalar *array_solver_x;

    PetscCall(VecGetOwnershipRange(mysolver->solver_x, &x_row_start, &x_row_end));
    PetscCall(VecGetLocalSize(mysolver->solver_x, &n_solver_x));
    PetscCall(VecGetArrayRead(mysolver->solver_x, &array_solver_x));

    FILE *fp = NULL;
    char path_sol[PETSC_MAX_PATH_LEN];
    sprintf(path_sol, "./answer_file/answer_x_%d.txt", myrank);

    if ((fp = fopen(path_sol, "wb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - \'%s\'\n", path_sol);
        exit(EXIT_FAILURE);
    }

    fprintf(fp, "%d\n", n_solver_x);
    for (int index = x_row_start; index < x_row_end; ++index)
    {
        fprintf(fp, "%d\t%021.16le\n", index, array_solver_x[index - x_row_start]);
    }

    fclose(fp);
}

void SolverPetscDestroy(MySolver *mysolver)
{
    PetscCall(KSPDestroy(&(mysolver->ksp)));
    PetscCall(MatDestroy(&(mysolver->solver_a)));
    PetscCall(VecDestroy(&(mysolver->solver_b)));
    PetscCall(VecDestroy(&(mysolver->solver_r)));
    PetscCall(VecDestroy(&(mysolver->solver_x)));
}
