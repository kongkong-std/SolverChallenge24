#include "main.h"

int main(int argc, char **argv)
{
    char *path_mat = NULL, *path_rhs = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-mat", argv[index]))
        {
            path_mat = argv[index + 1];
        }
        if (strstr("-rhs", argv[index]))
        {
            path_rhs = argv[index + 1];
        }
    }

    Matrix mat_a;
    Vector rhs_b;
    MatrixProcessSize(path_mat, &mat_a);
    VectorProcessSize(path_rhs, &rhs_b);

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    Mat solver_mat;
    Vec solver_rhs;
    int n_size = rhs_b.n;

    PetscCall(MatCreate(PETSC_COMM_WORLD, &solver_mat));
    PetscCall(MatSetSizes(solver_mat, PETSC_DECIDE, PETSC_DECIDE, n_size, n_size));
    PetscCall(MatSetType(solver_mat, MATAIJ));
    PetscCall(MatSetUp(solver_mat));

    int n_mat_loc_start = 0, n_mat_loc_end = 0;
    PetscCall(MatGetOwnershipRange(solver_mat, &n_mat_loc_start, &n_mat_loc_end));
    MatrixProcess(path_mat, &mat_a, n_mat_loc_start, n_mat_loc_end);

    for (int index = 0; index < mat_a.nnz; ++index)
    {
        int index_i = mat_a.row_idx[index];
        int index_j = mat_a.col_idx[index];
        PetscScalar val_tmp = mat_a.val_re[index] + mat_a.val_im[index] * PETSC_i;
        PetscCall(MatSetValue(solver_mat, index_i, index_j, val_tmp, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(solver_mat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(solver_mat, MAT_FINAL_ASSEMBLY));

#if 0
    PetscReal l1_norm_solver_mat = 0.;
    PetscCall(MatNorm(solver_mat, NORM_FROBENIUS, &l1_norm_solver_mat));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> L1 norm of matrix A = %021.16le\n", l1_norm_solver_mat));
#endif

    PetscCall(VecCreate(PETSC_COMM_WORLD, &solver_rhs));
    PetscCall(VecSetType(solver_rhs, VECMPI));
    PetscCall(VecSetSizes(solver_rhs, PETSC_DECIDE, n_size));

    int n_vec_loc_start = 0, n_vec_loc_end = 0;
    PetscCall(VecGetOwnershipRange(solver_rhs, &n_vec_loc_start, &n_vec_loc_end));
    VectorProcess(path_rhs, &rhs_b, n_vec_loc_start, n_vec_loc_end);

    for (int index = n_vec_loc_start; index < n_vec_loc_end; ++index)
    {
        PetscScalar val_tmp = rhs_b.val_re[index - n_vec_loc_start] + rhs_b.val_im[index - n_vec_loc_start] * PETSC_i;
        PetscCall(VecSetValues(solver_rhs, 1, &index, &val_tmp, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(solver_rhs));
    PetscCall(VecAssemblyEnd(solver_rhs));

#if 0
    PetscReal l2_norm_solver_rhs = 0.;
    PetscCall(VecNorm(solver_rhs, NORM_2, &l2_norm_solver_rhs));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, ">>>> L2 norm of vector b = %021.16le\n", l2_norm_solver_rhs));
#endif

#if 0
    PetscCall(MatView(solver_mat, PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(VecView(solver_rhs, PETSC_VIEWER_STDOUT_WORLD));
#endif

#if 1
    PetscViewer fd;
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "petsc_bin_mat", FILE_MODE_WRITE, &fd));
    PetscCall(MatView(solver_mat, fd));
    PetscCall(PetscViewerPopFormat(fd));
    PetscCall(PetscViewerDestroy(&fd));

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_WORLD, "petsc_bin_rhs", FILE_MODE_WRITE, &fd));
    PetscCall(VecView(solver_rhs, fd));
    PetscCall(PetscViewerPopFormat(fd));
    PetscCall(PetscViewerDestroy(&fd));
#endif

    // free memory
    free(mat_a.row_idx);
    free(mat_a.col_idx);
    free(mat_a.val_re);
    free(mat_a.val_im);
    free(rhs_b.val_re);
    free(rhs_b.val_im);

    PetscCall(MatDestroy(&solver_mat));
    PetscCall(VecDestroy(&solver_rhs));

    PetscFinalize();
    return 0;
}

/*
 * mpirun -np <np> ./app_agmg -mat ../input/real_mat.txt -rhs ../input/real_rhs.txt
 */