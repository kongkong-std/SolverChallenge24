#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct matrix
{
    int m, n, nnz;
    int *row_idx, *col_idx;
    double *val;
} Matrix;

typedef struct vector
{
    int n;
    double *val;
} Vector;

/*
 * matrix size data: m, n, nnz
 */
void MatrixProcessSize(const char *, Matrix *);

/*
 * vector size data: n
 */
void VectorProcessSize(const char *, Vector *);

void MatrixProcess(const char *, Matrix *, int, int);
void VectorProcess(const char *, Vector *, int, int);

#endif
