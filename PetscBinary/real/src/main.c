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
    MatrixProcess(path_mat, &mat_a);
    VectorProcess(path_rhs, &rhs_b);

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &argv, (char *)0, NULL));

    Mat solver_mat;
    Vec solver_rhs;
    int n_size = rhs_b.n;

    PetscCall(MatCreate(PETSC_COMM_SELF, &solver_mat));
    PetscCall(MatSetSizes(solver_mat, PETSC_DECIDE, PETSC_DECIDE, n_size, n_size));
    PetscCall(MatSetType(solver_mat, MATAIJ));
    PetscCall(MatSetUp(solver_mat));

    // 1-base to 0-base
    for (int index = 0; index < mat_a.nnz; ++index)
    {
        int index_i = mat_a.row_idx[index] - 1;
        int index_j = mat_a.col_idx[index] - 1;
        PetscScalar val_tmp = mat_a.val[index];
        PetscCall(MatSetValue(solver_mat, index_i, index_j, val_tmp, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(solver_mat, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(solver_mat, MAT_FINAL_ASSEMBLY));

    PetscCall(VecCreate(PETSC_COMM_SELF, &solver_rhs));
    PetscCall(VecSetType(solver_rhs, VECSEQ));
    PetscCall(VecSetSizes(solver_rhs, PETSC_DECIDE, n_size));

    for (int index = 0; index < n_size; ++index)
    {
        PetscScalar val_tmp = rhs_b.val[index];
        PetscCall(VecSetValues(solver_rhs, 1, &index, &val_tmp, INSERT_VALUES));
    }

    PetscCall(VecAssemblyBegin(solver_rhs));
    PetscCall(VecAssemblyEnd(solver_rhs));

    // petsc binary file
    PetscViewer fd;
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, "petsc_bin_mat", FILE_MODE_WRITE, &fd));
    PetscCall(MatView(solver_mat, fd));
    PetscCall(PetscViewerPopFormat(fd));
    PetscCall(PetscViewerDestroy(&fd));

    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, "petsc_bin_rhs", FILE_MODE_WRITE, &fd));
    PetscCall(VecView(solver_rhs, fd));
    PetscCall(PetscViewerPopFormat(fd));
    PetscCall(PetscViewerDestroy(&fd));

    // free memory
    free(mat_a.row_idx);
    free(mat_a.col_idx);
    free(mat_a.val);
    free(rhs_b.val);

    PetscCall(MatDestroy(&solver_mat));
    PetscCall(VecDestroy(&solver_rhs));

    PetscCall(PetscFinalize());
    return 0;
}

// command
/*
 * ./app_petsc_bin -mat <path/to/mat> -rhs <path/to/vec>
 */