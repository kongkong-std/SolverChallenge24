#include "utlise_long.h"

void spmv_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *x, long double *y)
{
    for (int i = 0; i < n; i++)
    {
        y[i] = 0.0;
        long double c = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++)
        {
            long double num = val[j] * x[col_idx[j]];
            long double z = num - c;
            long double t = y[i] + z;
            c = (t - y[i]) - z;
            y[i] = t;
        }
    }
}

// The complex number multiplication
void mul_ld(long double *v1, long double *v1i, long double *v2, long double *v2i, long double *v3, long double *v3i)
{
    *v3 = *v1 * *v2 - (*v1i * *v2i);
    *v3i = *v1i * *v2 + *v1 * *v2i;
}

// The complex number multiplication
void conj_mul_ld(long double *v1, long double *v1i, long double *v2, long double *v2i, long double *v3, long double *v3i)
{
    *v3 = *v1 * *v2 + (*v1i * *v2i);
    *v3i = *v1 * *v2i - *v1i * *v2;
}

// The complex number addition ,sum use kekan sum
void add_ld(long double *v1, long double *v1i, long double *sum, long double *sumi, long double *c, long double *ci)
{
    long double num = *v1;
    long double z = num - *c;
    long double t = *sum + z;
    *c = (t - *sum) - z;
    *sum = t;

    long double numi = *v1i;
    long double zi = numi - *ci;
    long double ti = *sumi + zi;
    *ci = (ti - *sumi) - zi;
    *sumi = ti;
}

// Multiply a csr matrix with two vector x and xi, and get the resulting vector y and yi
void spmv_complex_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *vali, long double *x, long double *xi, long double *y,
                     long double *yi)
{
    long double tmp[2];
    long double c[2];
    for (int i = 0; i < n; i++)
    {
        y[i] = yi[i] = 0.0;
        c[0] = c[1] = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++)
        {
            mul_ld(val + j, vali + j, x + col_idx[j], xi + col_idx[j], tmp, tmp + 1);
            add_ld(tmp, tmp + 1, y + i, yi + i, c, c + 1);
        }
    }
}

// Calculate the mode of complex number
void complex_modulus_squared_ld(long double *a, long double *ai, long double *l) { *l = *a * *a + (*ai * *ai); }

// Calculate the conjugate of complex number
void complex_conjugate_ld(long double *a, long double *ai, long double *b, long double *bi)
{
    *b = *a;
    *bi = -*ai;
}

// Calculate the division of complex number
void complex_division_ld(long double *a, long double *ai, long double *b, long double *bi, long double *c, long double *ci)
{
    long double tmp = 0.0;
    long double *l;
    l = &tmp;
    complex_modulus_squared_ld(b, bi, l);
    *c = *a * *b + (*ai * *bi);
    *ci = *ai * *b - (*a * *bi);
    *c = (*c) / (*l);
    *ci = (*ci) / (*l);
}

// Calculate the 2-norm of a vector ,sum use kekan sum
long double vec2norm_ld(long double *x, int n)
{
    long double sum = 0.0;
    long double c = 0.0;
    for (int i = 0; i < n; i++)
    {
        long double num = x[i] * x[i];
        long double z = num - c;
        long double t = sum + z;
        c = (t - sum) - z;
        sum = t;
    }

    return sqrtl(sum);
}

// Calculate the 2-norm of a complex vector ,sum use kekan sum
void vec2norm_complex_ld(long double *x, long double *xi, long double *err, long double *erri, int n)
{
    long double sum[2] = {0, 0}, tmp[2];
    long double c[2];
    c[0] = c[1] = 0.0;

    for (int i = 0; i < n; i++)
    {
        conj_mul_ld(x + i, xi + i, x + i, xi + i, tmp, tmp + 1);
        add_ld(tmp, tmp + 1, sum, sum + 1, c, c + 1);
    }

    *err = sqrtl(sum[0]);
    *erri = sqrtl(sum[1]);
}

long double max_check_ld(long double *x, int n)
{
    long double max = DBL_MIN;
    for (int i = 0; i < n; i++)
    {
        long double x_fabs = fabsl(x[i]);
        max = max > x_fabs ? max : x_fabs;
    }
    return max;
}

long double max_check_complex_ld(long double *x, long double *xi, int n)
{
    long double max = DBL_MIN;
    for (int i = 0; i < n; i++)
    {
        long double x_fabs = x[i] * x[i] + xi[i] * xi[i];
        max = max > x_fabs ? max : x_fabs;
    }
    return sqrtl(max);
}

// precision check using long double type
// validate the x
// answer1 = sqrtl(|| A*x - b ||)
// answer2 = || A*x - b || MAX
// answer3 = sqrtl(|| A*x - b ||/|| b ||)
// answer4 = MAX { |b - Ax|_i / |b_i| }
void check_correctness_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *x, long double *b)
{
    long double *b_new = (long double *)malloc(sizeof(long double) * n);
    long double *check_b = (long double *)malloc(sizeof(long double) * n);
    // long double* r_b     = (long double*)malloc(sizeof(long double) * n);
    spmv_ld(n, row_ptr, col_idx, val, x, b_new);
    for (int i = 0; i < n; i++)
    {
        check_b[i] = b_new[i] - b[i];
        // r_b[i]     = fabsl(check_b[i]) / MAX(fabsl(b[i]), 1e-20);
    }

    long double answer1 = vec2norm_ld(check_b, n);
    long double answer2 = max_check_ld(check_b, n);
    long double answer3 = answer1 / vec2norm_ld(b, n);
    // long double answer4 = max_check_ld(r_b, n);

    fprintf(stdout, "LD-Check || b - Ax || 2             =  %12.6Le\n", answer1);
    fprintf(stdout, "LD-Check || b - Ax || MAX           =  %12.6Le\n", answer2);
    fprintf(stdout, "LD-Check || b - Ax || 2 / || b || 2 =  %12.6Le\n", answer3);
    // fprintf(stdout, "LD-Check MAX { |b - Ax|_i / |b_i| } =  %12.6Le\n", answer4);

    free(b_new);
    free(check_b);
    // free(r_b);
}

void check_correctness_ld_d2ld(int n, int *row_ptr, int *col_idx, double *val, double *x, double *b)
{
    //! Step 1: data type transformation: double -> long double
    int nnz = row_ptr[n] - row_ptr[0];
    // printf("nnz = %d\n", nnz);
    long double *val_ld = (long double *)malloc(sizeof(long double) * nnz);
    long double *x_ld = (long double *)malloc(sizeof(long double) * n);
    long double *b_ld = (long double *)malloc(sizeof(long double) * n);

    for (int i = 0; i < nnz; i++)
        val_ld[i] = (long double)val[i];

    for (int i = 0; i < n; i++)
    {
        x_ld[i] = (long double)x[i];
        b_ld[i] = (long double)b[i];
    }

    //! Step 2: Check
    long double *b_new = (long double *)malloc(sizeof(long double) * n);
    long double *check_b = (long double *)malloc(sizeof(long double) * n);
    // long double* r_b     = (long double*)malloc(sizeof(long double) * n);
    spmv_ld(n, row_ptr, col_idx, val_ld, x_ld, b_new);
    for (int i = 0; i < n; i++)
    {
        check_b[i] = b_new[i] - b_ld[i];
        // r_b[i]     = fabsl(check_b[i]) / MAX(fabsl(b_ld[i]), 1e-20);
    }

    long double answer1 = vec2norm_ld(check_b, n);
    long double answer2 = max_check_ld(check_b, n);
    long double answer3 = answer1 / vec2norm_ld(b_ld, n);
    // long double answer4 = max_check_ld(r_b, n);

    fprintf(stdout, "LD-Check || b - Ax || 2             =  %12.6Le\n", answer1);
    fprintf(stdout, "LD-Check || b - Ax || MAX           =  %12.6Le\n", answer2);
    fprintf(stdout, "LD-Check || b - Ax || 2 / || b || 2 =  %12.6Le\n", answer3);
    // fprintf(stdout, "LD-Check MAX { |b - Ax|_i / |b_i| } =  %12.6Le\n", answer4);

    //! Step 3: free memory
    free(val_ld);
    free(x_ld);
    free(b_ld);

    free(b_new);
    free(check_b);
    // free(r_b);
}

void check_correctness_complex_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *vali, long double *x, long double *xi,
                                  long double *b, long double *bi)
{
    long double *b_new = (long double *)malloc(sizeof(long double) * n);
    long double *b_new_i = (long double *)malloc(sizeof(long double) * n);
    long double *check_b = (long double *)malloc(sizeof(long double) * n);
    long double *check_b_i = (long double *)malloc(sizeof(long double) * n);

    long double *r_b = (long double *)malloc(sizeof(long double) * n);
    long double *r_b_i = (long double *)malloc(sizeof(long double) * n);
    // long double* rb_mode = (long double*)malloc(sizeof(long double) * n);
    // long double  tmp     = 0.0;
    // long double* l;
    // l = &tmp;

    spmv_complex_ld(n, row_ptr, col_idx, val, vali, x, xi, b_new, b_new_i);
    for (int i = 0; i < n; i++)
    {
        check_b[i] = b_new[i] - b[i];
        check_b_i[i] = b_new_i[i] - bi[i];
        // complex_division_ld(&check_b[i], &check_b_i[i], &b[i], &bi[i], r_b, r_b_i);
        // complex_modulus_squared_ld(r_b, r_b_i, l);
        // rb_mode[i] = sqrtl(*l);
    }

    long double err_check_b = 0.0;
    long double err_check_bi = 0.0;
    long double err_b = 0.0;
    long double err_bi = 0.0;

    vec2norm_complex_ld(check_b, check_b_i, &err_check_b, &err_check_bi, n);
    vec2norm_complex_ld(b, bi, &err_b, &err_bi, n);
    long double max_answer = max_check_complex_ld(check_b, check_b_i, n);
    // long double max_answer2 = max_check_ld(rb_mode, n);

    fprintf(stdout, "LD-Check complex || b - Ax || 2             =  %12.6Le \n", err_check_b);
    fprintf(stdout, "LD-Check complex || b - Ax || MAX           =  %12.6Le \n", max_answer);
    fprintf(stdout, "LD-Check complex || b - Ax || 2 / || b || 2 =  %12.6Le \n", err_check_b / err_b);
    // fprintf(stdout, "LD-Check complex MAX { |b - Ax|_i / |b_i| } =  %12.6Le \n", max_answer2);

    free(b_new);
    free(b_new_i);
    free(check_b);
    free(check_b_i);
    free(r_b);
    free(r_b_i);
    // free(rb_mode);
}

void check_correctness_complex_ld_d2ld(int n, int *row_ptr, int *col_idx, double *val, double *vali, double *x, double *xi, double *b, double *bi)
{
    //! Step 1: data type transformation: double -> long double
    int nnz = row_ptr[n] - row_ptr[0];
    // printf("nnz = %d\n", nnz);
    long double *val_ld = (long double *)malloc(sizeof(long double) * nnz);
    long double *vali_ld = (long double *)malloc(sizeof(long double) * nnz);
    long double *x_ld = (long double *)malloc(sizeof(long double) * n);
    long double *xi_ld = (long double *)malloc(sizeof(long double) * n);
    long double *b_ld = (long double *)malloc(sizeof(long double) * n);
    long double *bi_ld = (long double *)malloc(sizeof(long double) * n);

    for (int i = 0; i < nnz; i++)
    {
        val_ld[i] = (long double)val[i];
        vali_ld[i] = (long double)vali[i];
    }

    for (int i = 0; i < n; i++)
    {
        x_ld[i] = (long double)x[i];
        xi_ld[i] = (long double)xi[i];
        b_ld[i] = (long double)b[i];
        bi_ld[i] = (long double)bi[i];
    }

    //! Step 2: Check
    long double *b_new = (long double *)malloc(sizeof(long double) * n);
    long double *b_new_i = (long double *)malloc(sizeof(long double) * n);
    long double *check_b = (long double *)malloc(sizeof(long double) * n);
    long double *check_b_i = (long double *)malloc(sizeof(long double) * n);

    long double *r_b = (long double *)malloc(sizeof(long double) * n);
    long double *r_b_i = (long double *)malloc(sizeof(long double) * n);
    // long double* rb_mode = (long double*)malloc(sizeof(long double) * n);
    // long double  tmp     = 0.0;
    // long double* l;
    // l = &tmp;

    spmv_complex_ld(n, row_ptr, col_idx, val_ld, vali_ld, x_ld, xi_ld, b_new, b_new_i);
    for (int i = 0; i < n; i++)
    {
        check_b[i] = b_new[i] - b_ld[i];
        check_b_i[i] = b_new_i[i] - bi_ld[i];
        // complex_division_ld(&check_b[i], &check_b_i[i], &b_ld[i], &bi_ld[i], r_b, r_b_i);
        // complex_modulus_squared_ld(r_b, r_b_i, l);
        // rb_mode[i] = sqrtl(*l);
    }

    long double err_check_b = 0.0;
    long double err_check_bi = 0.0;
    long double err_b = 0.0;
    long double err_bi = 0.0;

    vec2norm_complex_ld(check_b, check_b_i, &err_check_b, &err_check_bi, n);
    vec2norm_complex_ld(b_ld, bi_ld, &err_b, &err_bi, n);
    long double max_answer = max_check_complex_ld(check_b, check_b_i, n);
    // long double max_answer2 = max_check_ld(rb_mode, n);

    fprintf(stdout, "LD-Check complex || b - Ax || 2             =  %12.6Le \n", err_check_b);
    fprintf(stdout, "LD-Check complex || b - Ax || MAX           =  %12.6Le \n", max_answer);
    fprintf(stdout, "LD-Check complex || b - Ax || 2 / || b || 2 =  %12.6Le \n", err_check_b / err_b);
    // fprintf(stdout, "LD-Check complex MAX { |b - Ax|_i / |b_i| } =  %12.6Le \n", max_answer2);

    //! Step 3: free memory
    free(val_ld);
    free(vali_ld);
    free(x_ld);
    free(xi_ld);
    free(b_ld);
    free(bi_ld);

    free(b_new);
    free(b_new_i);
    free(check_b);
    free(check_b_i);
    free(r_b);
    free(r_b_i);
    // free(rb_mode);
}

// store x to a file
void store_x_ld(int n, long double *x, char *filename)
{
    FILE *p = fopen(filename, "w");
    fprintf(p, "%d\n", n);
    for (int i = 0; i < n; i++)
        fprintf(p, "%021.16Le\n", x[i]);
    fclose(p);
}

// load right-hand side vector b
void load_vector_ld(int n, long double *b, char *filename)
{
    FILE *p = fopen(filename, "r");
    int n_right;
    // MM_typecode matcode;

    char line[MM_MAX_LINE_LENGTH];
    fgets(line, MM_MAX_LINE_LENGTH, p);
    int r = 0;

    if (line[0] == '%')
    {
        r = fscanf(p, "%d", &n_right);
    }
    else
    {
        n_right = atoi(line);
    }

    if (n_right != n)
    {
        fclose(p);
        fprintf(stdout, "Invalid size of vector.\n");
        return;
    }
    for (int i = 0; i < n_right; i++)
        r = fscanf(p, "%Lf", &b[i]);
    fclose(p);
}
