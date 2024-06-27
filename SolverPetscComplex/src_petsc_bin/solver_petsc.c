#include "mysolver.h"
#include "solver_pcshell.h"

#ifdef CHALLENGE_06
void SolverRealPartMatrix(const Mat *mat, Mat *mat_re)
{
    PetscInt row_loc = 0, col_loc = 0;
    PetscCall(MatGetLocalSize(*mat, &row_loc, &col_loc));

    PetscBool done = PETSC_TRUE;
    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    PetscInt loc_row_start = 0, loc_row_end = 0;
    PetscCall(MatGetRowIJ(*mat, 0, PETSC_FALSE, PETSC_TRUE, NULL, &csr_ia, &csr_ja, &done));
    PetscCall(MatGetOwnershipRange(*mat, &loc_row_start, &loc_row_end));

    PetscCall(MatCreate(PETSC_COMM_WORLD, mat_re));
    PetscCall(MatSetSizes(*mat_re, row_loc, col_loc, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetType(*mat_re, MATAIJ));
    PetscCall(MatSetUp(*mat_re));

    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        int index_start = csr_ia[index - loc_row_start];
        int index_end = csr_ia[index - loc_row_start + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(*mat, index, csr_ja[index_j], &val_tmp));

            PetscScalar val_tmp_re = PetscRealPart(val_tmp);
            PetscCall(MatSetValues(*mat_re, 1, &index, 1, csr_ja + index_j, &val_tmp_re, INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(*mat_re, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*mat_re, MAT_FINAL_ASSEMBLY));
}

void SolverImaginaryPartMatrix(const Mat *mat, Mat *mat_im)
{
    PetscInt row_loc = 0, col_loc = 0;
    PetscCall(MatGetLocalSize(*mat, &row_loc, &col_loc));

    PetscBool done = PETSC_TRUE;
    const PetscInt *csr_ia = NULL;
    const PetscInt *csr_ja = NULL;
    PetscInt loc_row_start = 0, loc_row_end = 0;
    PetscCall(MatGetRowIJ(*mat, 0, PETSC_FALSE, PETSC_TRUE, NULL, &csr_ia, &csr_ja, &done));
    PetscCall(MatGetOwnershipRange(*mat, &loc_row_start, &loc_row_end));

    PetscCall(MatCreate(PETSC_COMM_WORLD, mat_im));
    PetscCall(MatSetSizes(*mat_im, row_loc, col_loc, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetType(*mat_im, MATAIJ));
    PetscCall(MatSetUp(*mat_im));

    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        int index_start = csr_ia[index - loc_row_start];
        int index_end = csr_ia[index - loc_row_start + 1];
        for (int index_j = index_start; index_j < index_end; ++index_j)
        {
            PetscScalar val_tmp;
            PetscCall(MatGetValue(*mat, index, csr_ja[index_j], &val_tmp));

            PetscScalar val_tmp_im = PetscImaginaryPart(val_tmp);
            PetscCall(MatSetValues(*mat_im, 1, &index, 1, csr_ja + index_j, &val_tmp_im, INSERT_VALUES));
        }
    }

    PetscCall(MatAssemblyBegin(*mat_im, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(*mat_im, MAT_FINAL_ASSEMBLY));
}

void SoverRealPartVector(const Vec *vec, Vec *vec_re)
{
    PetscInt row_loc = 0;
    PetscCall(VecGetLocalSize(*vec, &row_loc));

    PetscInt loc_row_start = 0, loc_row_end = 0;
    PetscCall(VecGetOwnershipRange(*vec, &loc_row_start, &loc_row_end));

    PetscCall(VecCreate(PETSC_COMM_WORLD, vec_re));
    PetscCall(VecSetSizes(*vec_re, row_loc, PETSC_DETERMINE));
    PetscCall(VecSetUp(*vec_re));

    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        PetscScalar val_tmp;
        PetscCall(VecGetValues(*vec, 1, &index, &val_tmp));

        PetscScalar val_tmp_re = PetscRealPart(val_tmp);
        PetscCall(VecSetValues(*vec_re, 1, &index, &val_tmp_re, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(*vec_re));
    PetscCall(VecAssemblyEnd(*vec_re));
}

void SolverImaginaryPartVector(const Vec *vec, Vec *vec_im)
{
    PetscInt row_loc = 0;
    PetscCall(VecGetLocalSize(*vec, &row_loc));

    PetscInt loc_row_start = 0, loc_row_end = 0;
    PetscCall(VecGetOwnershipRange(*vec, &loc_row_start, &loc_row_end));

    PetscCall(VecCreate(PETSC_COMM_WORLD, vec_im));
    PetscCall(VecSetSizes(*vec_im, row_loc, PETSC_DETERMINE));
    PetscCall(VecSetUp(*vec_im));

    for (int index = loc_row_start; index < loc_row_end; ++index)
    {
        PetscScalar val_tmp;
        PetscCall(VecGetValues(*vec, 1, &index, &val_tmp));

        PetscScalar val_tmp_im = PetscImaginaryPart(val_tmp);
        PetscCall(VecSetValues(*vec_im, 1, &index, &val_tmp_im, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(*vec_im));
    PetscCall(VecAssemblyEnd(*vec_im));
}
#endif // CHALLENGE_06

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
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "PETSc Binary Matrix file: %s\n", path_mat));
    }

    PetscCall(PetscOptionsGetString(NULL, NULL, "-file_rhs", path_rhs, sizeof(path_rhs), &path_flag));
    if (path_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "PETSc Binary RHS file: %s\n", path_rhs));
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

#ifdef CHALLENGE_06
    SolverRealPartMatrix(&(mysolver->solver_a), &(mysolver->solver_a_re));
    SolverImaginaryPartMatrix(&(mysolver->solver_a), &(mysolver->solver_a_im));
    SolverRealPartVector(&(mysolver->solver_b), &(mysolver->solver_b_re));
    SolverImaginaryPartVector(&(mysolver->solver_b), &(mysolver->solver_b_im));
    SolverRealPartVector(&(mysolver->solver_x), &(mysolver->solver_x_re));
    SolverImaginaryPartVector(&(mysolver->solver_x), &(mysolver->solver_x_im));
    SolverRealPartVector(&(mysolver->solver_r), &(mysolver->solver_r_re));
    SolverImaginaryPartVector(&(mysolver->solver_r), &(mysolver->solver_r_im));

    PetscCall(MatDuplicate(mysolver->solver_a_im, MAT_DO_NOT_COPY_VALUES, &(mysolver->solver_a_im_oppo)));
    PetscCall(MatAXPY(mysolver->solver_a_im_oppo, -1., mysolver->solver_a_im, SAME_NONZERO_PATTERN));

    Mat mat_array[4] = {mysolver->solver_a_re, mysolver->solver_a_im_oppo, mysolver->solver_a_im, mysolver->solver_a_re};
    PetscCall(MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, mat_array, &(mysolver->solver_block_a)));

    Vec rhs_vec_array[2] = {mysolver->solver_b_re, mysolver->solver_b_im};
    Vec sol_vec_array[2] = {mysolver->solver_x_re, mysolver->solver_x_im};
    Vec res_vec_array[2] = {mysolver->solver_r_re, mysolver->solver_r_im};

    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, rhs_vec_array, &(mysolver->solver_block_b)));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, sol_vec_array, &(mysolver->solver_block_x)));
    PetscCall(VecCreateNest(PETSC_COMM_WORLD, 2, NULL, res_vec_array, &(mysolver->solver_block_r)));
#endif // CHALLENGE_06
}

void SolverPetscPreprocess(int argc, char **argv, MySolver *mysolver)
{
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &(mysolver->ksp)));
#ifdef CHALLENGE_06
    PetscReal shift_pc_re = 0., shift_pc_im = 0.;
    PetscBool shift_flag;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-shift_pc_re", &shift_pc_re, &shift_flag));
    if (shift_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "value for shift of preconditioner: %021.16le\n", shift_pc_re));
    }
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-shift_pc_im", &shift_pc_im, &shift_flag));
    if (shift_flag)
    {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "value for shift of preconditioner: %021.16le\n", shift_pc_im));
    }

    // shift pc
    /*
     * solver_a = solver_a_re + solver_a_im * i
     * solver_pc = solver_pc_re + solver_pc_im * i
     *     solver_pc_re = solver_a_re + shift_re * solver_a_im
     *     solver_pc_im = (1 + shift_im) * solver_a_im
     *
     * solver_pc_block = [solver_pc_re    -solver_pc_im]
     *                   [solver_pc_im     solver_pc_re]
     */
    Mat solver_pc_re, solver_pc_im, solver_pc_im_oppo;
    PetscCall(MatDuplicate(mysolver->solver_a_re, MAT_DO_NOT_COPY_VALUES, &solver_pc_re));
    PetscCall(MatDuplicate(mysolver->solver_a_im, MAT_DO_NOT_COPY_VALUES, &solver_pc_im));
    PetscCall(MatDuplicate(mysolver->solver_a_im, MAT_DO_NOT_COPY_VALUES, &solver_pc_im_oppo));

    PetscCall(MatAXPY(solver_pc_re, shift_pc_re, mysolver->solver_a_im, SAME_NONZERO_PATTERN));
    PetscCall(MatAXPY(solver_pc_re, 1., mysolver->solver_a_re, SAME_NONZERO_PATTERN));
    PetscReal shift_pc_im_tmp = 1. + shift_pc_im;
    PetscCall(MatAXPY(solver_pc_im, shift_pc_im_tmp, mysolver->solver_a_im, SAME_NONZERO_PATTERN));
    PetscCall(MatAXPY(solver_pc_im_oppo, -1., solver_pc_im, SAME_NONZERO_PATTERN));

    Mat pc_array[4] = {solver_pc_re, solver_pc_im_oppo, solver_pc_im, solver_pc_re};
    Mat solver_block_pc;
    PetscCall(MatCreateNest(PETSC_COMM_WORLD, 2, NULL, 2, NULL, pc_array, &solver_block_pc));

    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, solver_block_pc));
#elif
    PetscCall(KSPSetOperators(mysolver->ksp, mysolver->solver_a, mysolver->solver_a));
#endif
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
#ifdef CHALLENGE_06
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_block_b, mysolver->solver_block_x));
#elif
    PetscCall(KSPSolve(mysolver->ksp, mysolver->solver_b, mysolver->solver_x));
#endif
}

void SolverPetscResidualCheck(int argc, char **argv, MySolver *mysolver)
{
    PetscReal b_norm_2 = 0.;
    PetscReal r_norm_2 = 0.;

#ifdef CHALLENGE_06
    PetscCall(VecNorm(mysolver->solver_block_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_block_a, mysolver->solver_block_x, mysolver->solver_block_r));
    PetscCall(VecAXPY(mysolver->solver_block_r, -1., mysolver->solver_block_b));
    PetscCall(VecNorm(mysolver->solver_block_r, NORM_2, &r_norm_2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "            || b ||_2 = %021.16le\n", b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "|| r ||_2 / || b ||_2 = %021.16le\n", r_norm_2 / b_norm_2));
#elif
    PetscCall(VecNorm(mysolver->solver_b, NORM_2, &b_norm_2));

    PetscCall(MatMult(mysolver->solver_a, mysolver->solver_x, mysolver->solver_r));
    PetscCall(VecAXPY(mysolver->solver_r, -1., mysolver->solver_b));
    PetscCall(VecNorm(mysolver->solver_r, NORM_2, &r_norm_2));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "            || b ||_2 = %021.16le\n", b_norm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "|| r ||_2 / || b ||_2 = %021.16le\n", r_norm_2 / b_norm_2));
#endif

    PetscCall(KSPDestroy(&(mysolver->ksp)));
    PetscCall(MatDestroy(&(mysolver->solver_a)));
    PetscCall(VecDestroy(&(mysolver->solver_b)));
    PetscCall(VecDestroy(&(mysolver->solver_r)));
    PetscCall(VecDestroy(&(mysolver->solver_x)));
#ifdef CHALLENGE_06
    PetscCall(MatDestroy(&(mysolver->solver_a_re)));
    PetscCall(MatDestroy(&(mysolver->solver_a_im)));
    PetscCall(MatDestroy(&(mysolver->solver_a_im_oppo)));
    PetscCall(VecDestroy(&(mysolver->solver_b_re)));
    PetscCall(VecDestroy(&(mysolver->solver_b_im)));
    PetscCall(VecDestroy(&(mysolver->solver_x_re)));
    PetscCall(VecDestroy(&(mysolver->solver_x_im)));
    PetscCall(VecDestroy(&(mysolver->solver_r_re)));
    PetscCall(VecDestroy(&(mysolver->solver_r_im)));
#endif
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
