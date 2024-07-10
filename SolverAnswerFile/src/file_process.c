#include "main.h"

void AnswerXFileProcess(const char *path, Vector *vec)
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

    if ((vec->row_idx = (int *)malloc(n * sizeof(int))) == NULL ||
        (vec->val_re = (double *)malloc(n * sizeof(double))) == NULL ||
        (vec->val_im = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'answer file\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        fscanf(fp, "%d%lf%lf", vec->row_idx + index,
               vec->val_re + index,
               vec->val_im + index);
    }

    fclose(fp);
}