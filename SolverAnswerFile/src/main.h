#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUFFER_SIZE 4096

// #define REAL_VECTOR_
// #define COMPLEX_VECTOR_

#ifdef COMPLEX_VECTOR_
typedef struct vector
{
    /* data */
    int n;
    int *row_idx;
    double *val_re, *val_im;
} Vector;
#elif defined REAL_VECTOR_
typedef struct vector
{
    /* data */
    int n;
    int *row_idx;
    double *val;
} Vector;

#endif

void AnswerXFileProcess(const char *, Vector *);

#endif // MAIN_H_