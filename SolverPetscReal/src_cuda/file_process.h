#ifndef FILE_PROCESS_H_
#define FILE_PROCESS_H_

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



void MatrixProcess(const char *, Matrix *);
void VectorProcess(const char *, Vector *);

#endif // FILE_PROCESS_H_
