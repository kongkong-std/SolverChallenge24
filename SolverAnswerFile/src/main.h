#ifndef MAIN_H_
#define MAIN_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_BUFFER_SIZE 4096

typedef struct vector
{
    /* data */
    int n;
    int *row_idx;
    double *val_re, *val_im;
} Vector;

void AnswerXFileProcess(const char *, Vector *);

#endif    // MAIN_H_