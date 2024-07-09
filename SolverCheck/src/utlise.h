// Thanks to Yi Zong for providing support!
// check memory using function
#ifndef UTLISE_H_
#define UTLISE_H_

// user define by shb & milkyway
#define MM_MAX_LINE_LENGTH __INT32_MAX__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <math.h>
#include <sys/time.h>
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

void mem_usage();

// check answer function
void check_correctness(int n, int *row_ptr, int *col_idx, double *val, double *x, double *b);

void check_correctness_complex(int n, int *row_ptr, int *col_idx, double *val, double *val_v, double *x, double *xi, double *b, double *bi);

// store vector function
void store_x(int n, double *x, char *filename);

void print_help();

// Multiply a csr matrix with a vector x, and get the resulting vector y ,sum use kekan sum
void spmv(int n, int *row_ptr, int *col_idx, double *val, double *x, double *y);

// The complex number multiplication
void mul(double *v1, double *v1i, double *v2, double *v2i, double *v3, double *v3i);

// The complex number multiplication
void conj_mul(double *v1, double *v1i, double *v2, double *v2i, double *v3, double *v3i);
// The complex number addition ,sum use kekan sum
void add(double *v1, double *v1i, double *sum, double *sumi, double *c, double *ci);

// Multiply a csr matrix with two vector x and xi, and get the resulting vector y and yi
void spmv_complex(int n, int *row_ptr, int *col_idx, double *val,
                  double *vali, double *x, double *xi, double *y, double *yi);

// Calculate the mode of complex number
void complex_modulus_squared(double *a, double *ai, double *l);

// Calculate the conjugate of complex number
void complex_conjugate(double *a, double *ai, double *b, double *bi);

// Calculate the division of complex number
void complex_division(double *a, double *ai, double *b, double *bi, double *c, double *ci);

// Calculate the 2-norm of a vector ,sum use kekan sum
double vec2norm(double *x, int n);

// Calculate the 2-norm of a complex vector ,sum use kekan sum
void vec2norm_complex(double *x, double *xi, double *err, double *erri, int n);

double max_check(double *x, int n);

double max_check_complex(double *x, double *xi, int n);

// store x (complex type) to a file
void store_x_complex(int n, double *x, double *x_v, char *filename);

// return now time
double GetCurrentTime();

#endif