// Thanks to Haibing Sun and Li Zhao for providing support for long double precision
// detection. Here, long double is abbreviated as "ld"

#ifndef UTLISE_LONG_H_
#define UTLISE_LONG_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>

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
void check_correctness_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *x, long double *b);

void check_correctness_ld_d2ld(int n, int *row_ptr, int *col_idx, double *val, double *x, double *b);

void check_correctness_complex_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *vali, long double *x, long double *xi,
                                  long double *b, long double *bi);

void check_correctness_complex_ld_d2ld(int n, int *row_ptr, int *col_idx, double *val, double *vali, double *x, double *xi, double *b, double *bi);

// store vector function
void store_x_ld(int n, long double *x, char *filename);

// Multiply a csr matrix with a vector x, and get the resulting vector y ,sum use kekan
// sum
void spmv_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *x, long double *y);

// The complex number multiplication
void mul_ld(long double *v1, long double *v1i, long double *v2,
            long double *v2i, long double *v3, long double *v3i);

// The complex number multiplication
void conj_mul_ld(long double *v1, long double *v1i, long double *v2, long double *v2i,
                 long double *v3, long double *v3i);

// The complex number addition ,sum use kekan sum
void add_ld(long double *v1, long double *v1i, long double *sum, long double *sumi,
            long double *c, long double *ci);

// Multiply a csr matrix with two vector x and xi, and get the resulting vector y and yi
void spmv_complex_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *vali, long double *x, long double *xi, long double *y,
                     long double *yi);

// Calculate the mode of complex number
void complex_modulus_squared_ld(long double *a, long double *ai, long double *l);

// Calculate the conjugate of complex number
void complex_conjugate_ld(long double *a, long double *ai, long double *b, long double *bi);
// Calculate the division of complex number
void complex_division_ld(long double *a, long double *ai, long double *b, long double *bi,
                         long double *c, long double *ci);

// Calculate the 2-norm of a vector ,sum use kekan sum
long double vec2norm_ld(long double *x, int n);

// Calculate the 2-norm of a complex vector ,sum use kekan sum
void vec2norm_complex_ld(long double *x, long double *xi, long double *err, long double *erri, int n);
long double max_check_ld(long double *x, int n);

long double max_check_complex_ld(long double *x, long double *xi, int n);

// precision check using long double type
// validate the x
// answer1 = sqrtl(|| A*x - b ||)
// answer2 = || A*x - b || MAX
// answer3 = sqrtl(|| A*x - b ||/|| b ||)
// answer4 = MAX { |b - Ax|_i / |b_i| }
void check_correctness_ld(int n, int *row_ptr, int *col_idx, long double *val,
                          long double *x, long double *b);

void check_correctness_ld_d2ld(int n, int *row_ptr, int *col_idx, double *val, double *x, double *b);

void check_correctness_complex_ld(int n, int *row_ptr, int *col_idx, long double *val, long double *vali, long double *x, long double *xi,
                                  long double *b, long double *bi);

void check_correctness_complex_ld_d2ld(int n, int *row_ptr, int *col_idx,
                                       double *val, double *vali, double *x, double *xi, double *b, double *bi);

// store x to a file
void store_x_ld(int n, long double *x, char *filename);
#endif