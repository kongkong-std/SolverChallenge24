// Thanks to Yi Zong for providing support!
// check memory using function

// user define by shb & milkyway
#define MM_MAX_LINE_LENGTH  __INT32_MAX__
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <float.h>

/**
 * \brief Definition of max, min, abs
 */
#ifndef _MAX_MIN_ABS_
#define _MAX_MIN_ABS_
#define MAX(a, b) (((a) > (b)) ? (a) : (b)) ///< bigger one in a and b
#define MIN(a, b) (((a) < (b)) ? (a) : (b)) ///< smaller one in a and b
#define ABS(a) (((a) >= 0.0) ? (a) : -(a))  ///< absolute value of a
#endif

// check answer function
void check_correctness(int n, int* row_ptr, int* col_idx, double* val, double* x, double* b);

void check_correctness_complex(int n, int* row_ptr, int* col_idx, double* val, double* val_v, double* x, double* xi, double* b, double* bi);

// store vector function
void store_x(int n, double* x, char* filename);

// read vector function
void load_vector(int n, double* b, char* filename);

void mem_usage();
void mem_usage()
{
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    long long int mymem = r_usage.ru_maxrss;
#ifdef MPI_USE
    long long int total;
    MPI_Allreduce(&mymem, &total, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0)
    {
        fprintf(stdout, "Total memory usage:                 %.4f MB\n", (double)((double)(total) / (double)(1024)));
        fflush(stdout);
    }
#else
    fprintf(stdout, "Total memory usage:                 %.4f MB\n", (double)((double)(mymem) / (double)(1024)));
    fflush(stdout);
#endif
}

void print_help()
{
    fprintf(stdout, "please follow this input:\n./solver matrixname.mtx vectorb.rhs sol.rhs\n");
    fprintf(stdout, "e.g., ./solver matrixname.mtx vectorb.rhs 1 0 0\n");
    fprintf(stdout, "\nhelp:\n");
    fprintf(stdout, "matrixname.mtx: the filename of matrix A\n");
    fprintf(stdout, "vectorb.rhs:    the filename of right-hand side vector b\n");
    fprintf(stdout, "read_base:      0-base or 1-base, default 1-base\n");
    fprintf(stdout, "type:           type to output time, 0: end to end time, default 0\n");
    fprintf(stdout, "                                     1: solver time + solve time\n");
    fprintf(stdout, "                                     2: solve time\n");
    fprintf(stdout, "sys_type:       type of algebraic systems, 0: real, 1: complex; default 0\n");
}

// Multiply a csr matrix with a vector x, and get the resulting vector y ,sum use kekan sum
void spmv(int n, int* row_ptr, int* col_idx, double* val, double* x, double* y)
{
    for (int i = 0; i < n; i++) {
        y[i]     = 0.0;
        double c = 0.0;
        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            double num = val[j] * x[col_idx[j]];
            double z   = num - c;
            double t   = y[i] + z;
            c          = (t - y[i]) - z;
            y[i]       = t;
        }
    }
}

// The complex number multiplication
void mul(double* v1, double* v1i, double* v2, double* v2i, double* v3, double* v3i)
{
    *v3  = *v1 * *v2 - (*v1i * *v2i);
    *v3i = *v1i * *v2 + *v1 * *v2i;
}

// The complex number multiplication
void conj_mul(double* v1, double* v1i, double* v2, double* v2i, double* v3, double* v3i)
{
    *v3  = *v1 * *v2 + (*v1i * *v2i);
    *v3i = *v1 * *v2i - *v1i * *v2;
}

// The complex number addition ,sum use kekan sum
void add(double* v1, double* v1i, double* sum, double* sumi, double* c, double* ci)
{
    double num = *v1;
    double z   = num - *c;
    double t   = *sum + z;
    *c         = (t - *sum) - z;
    *sum       = t;

    double numi = *v1i;
    double zi   = numi - *ci;
    double ti   = *sumi + zi;
    *ci         = (ti - *sumi) - zi;
    *sumi       = ti;
}

// Multiply a csr matrix with two vector x and xi, and get the resulting vector y and yi
void spmv_complex(int n, int* row_ptr, int* col_idx, double* val, double* vali, double* x, double* xi, double* y, double* yi)
{
    double tmp[2];
    double c[2];
    for (int i = 0; i < n; i++) {
        y[i] = yi[i] = 0.0;
        c[0] = c[1] = 0.0;
        // printf(c + 1);

        for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            mul(val + j, vali + j, x + col_idx[j], xi + col_idx[j], tmp, tmp + 1);
            add(tmp, tmp + 1, y + i, yi + i, c, c + 1);
        }
    }
}

// Calculate the mode of complex number
void complex_modulus_squared(double* a, double* ai, double* l) { *l = *a * *a + (*ai * *ai); }

// Calculate the conjugate of complex number
void complex_conjugate(double* a, double* ai, double* b, double* bi)
{
    *b  = *a;
    *bi = -*ai;
}

// Calculate the division of complex number
void complex_division(double* a, double* ai, double* b, double* bi, double* c, double* ci)
{
    double  tmp = 0.0;
    double* l;
    l = &tmp;
    complex_modulus_squared(b, bi, l);
    *c  = *a * *b + (*ai * *bi);
    *ci = *ai * *b - (*a * *bi);
    *c  = (*c) / (*l);
    *ci = (*ci) / (*l);
}

// Calculate the 2-norm of a vector ,sum use kekan sum
double vec2norm(double* x, int n)
{
    double sum = 0.0;
    double c   = 0.0;
    for (int i = 0; i < n; i++) {
        double num = x[i] * x[i];
        double z   = num - c;
        double t   = sum + z;
        c          = (t - sum) - z;
        sum        = t;
    }

    return sqrt(sum);
}

// Calculate the 2-norm of a complex vector ,sum use kekan sum
void vec2norm_complex(double* x, double* xi, double* err, double* erri, int n)
{
    double sum[2] = {0, 0}, tmp[2];
    double c[2];
    c[0] = c[1] = 0.0;

    for (int i = 0; i < n; i++) {
        conj_mul(x + i, xi + i, x + i, xi + i, tmp, tmp + 1);
        add(tmp, tmp + 1, sum, sum + 1, c, c + 1);
    }

    *err  = sqrt(sum[0]);
    *erri = sqrt(sum[1]);
}

#include <float.h>
double max_check(double* x, int n)
{
    double max = DBL_MIN;
    for (int i = 0; i < n; i++) {
        double x_fabs = fabs(x[i]);
        max           = max > x_fabs ? max : x_fabs;
    }
    return max;
}

double max_check_complex(double* x, double* xi, int n)
{
    double max = DBL_MIN;
    for (int i = 0; i < n; i++) {
        double x_fabs = x[i] * x[i] + xi[i] * xi[i];
        max           = max > x_fabs ? max : x_fabs;
    }
    return sqrt(max);
}

// validate the x
// answer1 = sqrt(|| A*x - b ||)
// answer2 = || A*x - b || MAX
// answer3 = sqrt(|| A*x - b ||/|| b ||)
// precision check using double type
// void check_correctness(int n, int *row_ptr, int *col_idx, double *val, double *x, double *b)
void check_correctness(int n, int* row_ptr, int* col_idx, double* val, double* x, double* b)
{

    double* b_new   = (double*)malloc(sizeof(double) * n);
    double* check_b = (double*)malloc(sizeof(double) * n);
    // double* r_b     = (double*)malloc(sizeof(double) * n);

    spmv(n, row_ptr, col_idx, val, x, b_new);
    for (int i = 0; i < n; i++) {
        check_b[i] = b_new[i] - b[i];
        // r_b[i]     = fabs(check_b[i]) / MAX((b[i]), 1e-20);
    }

    double answer1 = vec2norm(check_b, n);
    double answer2 = max_check(check_b, n);
    double answer3 = answer1 / vec2norm(b, n);
    // double answer4 = max_check(r_b, n);
    fprintf(stdout, "Check || b - Ax || 2             =  %12.6e\n", answer1);
    fprintf(stdout, "Check || b - Ax || MAX           =  %12.6e\n", answer2);
    fprintf(stdout, "Check || b - Ax || 2 / || b || 2 =  %12.6e\n", answer3);
    // fprintf(stdout, "Check MAX { |b - Ax|_i / |b_i| } =  %12.6e\n", answer4);

    free(b_new);
    free(check_b);
    // free(r_b);
}

void check_correctness_complex(int n, int* row_ptr, int* col_idx, double* val, double* val_v, double* x, double* xi, double* b, double* bi)
{
    double* b_new     = (double*)malloc(sizeof(double) * n);
    double* b_new_i   = (double*)malloc(sizeof(double) * n);
    double* check_b   = (double*)malloc(sizeof(double) * n);
    double* check_b_i = (double*)malloc(sizeof(double) * n);

    double* r_b     = (double*)malloc(sizeof(double) * n);
    double* r_b_i   = (double*)malloc(sizeof(double) * n);
    // double* rb_mode = (double*)malloc(sizeof(double) * n);
    // double  tmp     = 0.0;
    // double* l;
    // l = &tmp;

    spmv_complex(n, row_ptr, col_idx, val, val_v, x, xi, b_new, b_new_i);
    for (int i = 0; i < n; i++) {
        check_b[i]   = b_new[i] - b[i];
        check_b_i[i] = b_new_i[i] - bi[i];
        // complex_division(&check_b[i], &check_b_i[i], &b[i], &bi[i], r_b, r_b_i);
        // complex_modulus_squared(r_b, r_b_i, l);
        // rb_mode[i] = sqrtl(*l);
    }

    double err_check_b  = 0.0;
    double err_check_bi = 0.0;
    double err_b        = 0.0;
    double err_bi       = 0.0;

    vec2norm_complex(check_b, check_b_i, &err_check_b, &err_check_bi, n);
    vec2norm_complex(b, bi, &err_b, &err_bi, n);
    double max_answer  = max_check_complex(check_b, check_b_i, n);
    // double max_answer2 = max_check(rb_mode, n);

    fprintf(stdout, "Check complex || b - Ax || 2             =  %12.6e \n", err_check_b);
    fprintf(stdout, "Check complex || b - Ax || MAX           =  %12.6e \n", max_answer);
    fprintf(stdout, "Check complex || b - Ax || 2 / || b || 2 =  %12.6e \n", err_check_b / err_b);
    // fprintf(stdout, "Check complex MAX { |b - Ax|_i / |b_i| } =  %12.6e \n", max_answer2);

    free(b_new);
    free(b_new_i);
    free(check_b);
    free(check_b_i);
    free(r_b);
    free(r_b_i);
    // free(rb_mode);
}

// store x to a file
void store_x(int n, double* x, char* filename)
{
    FILE* p = fopen(filename, "w");
    fprintf(p, "%d\n", n);
    for (int i = 0; i < n; i++) fprintf(p, "%24.20e\n", x[i]);
    fclose(p);
}

// store x (complex type) to a file
void store_x_complex(int n, double* x, double* x_v, char* filename)
{
    FILE* p = fopen(filename, "w");
    fprintf(p, "%d\n", n);
    int i;
    for (i = 0; i < n; i++) fprintf(p, "%lf %lf\n", x[i], x_v[i]);
    fclose(p);
}

// load right-hand side vector b
void load_vector(int n, double* b, char* filename)
{
    FILE* p = fopen(filename, "r");
    int   n_right;
    // MM_typecode matcode;

    int  r = 0;
    char line[MM_MAX_LINE_LENGTH];
    fgets(line, MM_MAX_LINE_LENGTH, p);

    if (line[0] == '%') {
        r = fscanf(p, "%d", &n_right);
    } else {
        n_right = atoi(line);
    }

    if (n_right != n) {
        fclose(p);
        fprintf(stdout, "Invalid size of vector.\n");
        return;
    }
    for (int i = 0; i < n_right; i++) r = fscanf(p, "%lf", &b[i]);
    fclose(p);
}

// load right-hand side vector b (complex type)
void load_b_complex(int n, double* b, double* b_v, char* filename)
{
    FILE* p = fopen(filename, "r");
    int   n_right;
    // MM_typecode matcode;

    int  r = 0;
    char line[MM_MAX_LINE_LENGTH];
    fgets(line, MM_MAX_LINE_LENGTH, p);

    if (line[0] == '%') {
        r = fscanf(p, "%d", &n_right);
    } else {
        n_right = atoi(line);
    }

    if (n_right != n) {
        fclose(p);
        fprintf(stdout, "Invalid size of vector.\n");
        return;
    }
    for (int i = 0; i < n_right; i++) r = fscanf(p, "%lf %lf", &b[i], &b_v[i]);
    fclose(p);
}

// return now time
double GetCurrentTime()
{
    struct timeval time;
    gettimeofday(&time, NULL);
    return time.tv_sec * 1000.0 + time.tv_usec / 1000.0;
}
