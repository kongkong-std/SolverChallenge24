#include "main.h"

int main(int argc, char **argv)
{
    char *path_mat = NULL, *path_rhs = NULL, *path_sol = NULL;
    char *dst_mat = NULL;
    int sys_type = 0; // 0, real system; 1, complex system
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-path_mat", argv[index]))
        {
            path_mat = argv[index + 1];
        }
        if (strstr("-path_rhs", argv[index]))
        {
            path_rhs = argv[index + 1];
        }
        if (strstr("-path_sol", argv[index]))
        {
            path_sol = argv[index + 1];
        }
        if (strstr("-sys_type", argv[index]))
        {
            sys_type = atoi(argv[index + 1]);
        }
        if (strstr("-dst_mat", argv[index]))
        {
            dst_mat = argv[index + 1];
        }
    }

    int m = 0, n = 0, nnz = 0;
    int *row_ptr = NULL, *col_idx = NULL;
    double *val = NULL, *val_im = NULL;
    double *b = NULL, *b_im = NULL;
    double *x = NULL, *x_im = NULL;

    // real system
    if (sys_type == 0)
    {
        RealCOO2CSRMatrixFileProcess(path_mat, &m, &n, &nnz, &row_ptr, &col_idx, &val);

        if (m != n)
        {
            fprintf(stderr, "Invalid matrix size.\n");
            exit(EXIT_FAILURE);
        }

        FILE *fp = NULL;
        if ((fp = fopen(dst_mat, "wb")) == NULL)
        {
            fprintf(stderr, "Cannot open file!\n");
            exit(EXIT_FAILURE);
        }

        fprintf(fp, "%d\n", m);
        for (int index = 0; index < m + 1; ++index)
        {
            fprintf(fp, "%d\n", row_ptr[index] + 1);
        }
        for (int index = 0; index < nnz; ++index)
        {
            fprintf(fp, "%d\n", col_idx[index] + 1);
        }
        for (int index = 0; index < nnz; ++index)
        {
            fprintf(fp, "%021.16le\n", val[index]);
        }

        fclose(fp);

#ifdef CHECK_ERROR_
        RealRHSFileProcess(path_rhs, &b);
        RealRHSFileProcess(path_sol, &x);

        // 1) using double precision
        check_correctness(n, row_ptr, col_idx, val, x, b);
        // 2) using long double precision
        check_correctness_ld_d2ld(n, row_ptr, col_idx, val, x, b);
        free(x);
        free(b);
#endif

        free(row_ptr);
        free(col_idx);
        free(val);
    }

    // complex system
    if (sys_type == 1)
    {
        ComplexCOO2CSRMatrixFileProcess(path_mat, &m, &n, &nnz, &row_ptr, &col_idx, &val, &val_im);

        if (m != n)
        {
            fprintf(stderr, "Invalid matrix size.\n");
            exit(EXIT_FAILURE);
        }

#ifdef CHECK_ERROR_
        ComplexRHSFileProcess(path_rhs, &b, &b_im);
        ComplexRHSFileProcess(path_sol, &x, &x_im);

        // 1) using double precision
        check_correctness_complex(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);
        // 2) using long double precision
        check_correctness_complex_ld_d2ld(n, row_ptr, col_idx, val, val_im, x, x_im, b, b_im);

        free(x);
        free(x_im);
        free(b);
        free(b_im);
#endif

        free(row_ptr);
        free(col_idx);
        free(val);
        free(val_im);
    }

    return 0;
}

// command
/*
 * -path_mat    <path/to/matrix>
 * -path_rhs    <path/to/rhs>
 * -path_sol    <path/to/sol>
 * -sys_type    <0/1>
 */
