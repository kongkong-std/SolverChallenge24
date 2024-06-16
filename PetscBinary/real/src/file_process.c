#include "main.h"
#define MAX_SIZE 512

void MatrixProcess(const char *path, Matrix *mat)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open matrix file - \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    char buffer[MAX_SIZE];

    // previous 2 lines ignored
    fgets(buffer, MAX_SIZE, fp);
    fgets(buffer, MAX_SIZE, fp);

    int m = 0, n = 0, nnz = 0;
    fscanf(fp, "%d%d%d", &m, &n, &nnz);
    mat->m = m;
    mat->n = n;
    mat->nnz = nnz;

    if ((mat->row_idx = (int *)malloc(nnz * sizeof(int))) == NULL ||
        (mat->col_idx = (int *)malloc(nnz * sizeof(int))) == NULL ||
        (mat->val = (double *)malloc(nnz * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, "%d%d%lf", mat->row_idx + index,
               mat->col_idx + index,
               mat->val + index);
    }

    fclose(fp);
}

void VectorProcess(const char *path, Vector *vec)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    int n = 0;
    fscanf(fp, "%d", &n);
    vec->n = n;

    if ((vec->val = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        exit(EXIT_FAILURE);
    }

    for(int index = 0; index < n; ++index)
    {
        fscanf(fp, "%d", vec->val + index);
    }

    fclose(fp);
}