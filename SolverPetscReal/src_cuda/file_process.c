#include "file_process.h"
#define MAX_SIZE 4096

void MatrixProcess(const char *path, Matrix *mat)
{
    int m = 0, n = 0, nnz = 0;

    FILE *fp = NULL;

    // scanning file to get local nnz
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix process \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    // skip first 2 lines
    char buffer[MAX_SIZE];
    fgets(buffer, MAX_SIZE, fp);
    fgets(buffer, MAX_SIZE, fp);

    fscanf(fp, "%d%d%d\n", &m, &n, &nnz);
    mat->m = m;
    mat->n = n;
    mat->nnz = nnz;

    // getting local row_idx, col_idx, val
    if ((mat->row_idx = (int *)malloc(nnz * sizeof(int))) == NULL ||
        (mat->col_idx = (int *)malloc(nnz * sizeof(int))) == NULL ||
        (mat->val = (double *)malloc(nnz * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'matrix data\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < nnz; ++index)
    {
        fscanf(fp, "%d%d%lf", mat->row_idx + index, mat->col_idx + index, mat->val + index);
    }

    fclose(fp);
}

void VectorProcess(const char *path, Vector *vec)
{
    int n = 0;

    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - vector process \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    fscanf(fp, "%d", &n);
    vec->n = n;

    if ((vec->val = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'vector data\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        fscanf(fp, "%lf", vec->val + index);
    }

    fclose(fp);
}
