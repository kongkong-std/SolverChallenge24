#ifndef MAIN_H_
#define MAIN_H_

#include <petscksp.h>

typedef struct matrix
{
    int m, n, nnz;
    int *row_idx, *col_idx;
    double *val;
} Matrix;

typedef struct vector
{
    int n;
    double * val;
} Vector;

void MatrixProcess(const char *, Matrix *);
void VectorProcess(const char *, Vector *);

#endif