#include "file_process.h"
#define MAX_SIZE 512

void MatrixProcessSize(const char *path, Matrix *mat)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix size \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    // skip first 2 lines
    char buffer[MAX_SIZE];
    fgets(buffer, MAX_SIZE, fp);
    fgets(buffer, MAX_SIZE, fp);

    fscanf(fp, "%d%d%d", &(mat->m), &(mat->n), &(mat->nnz));
    fclose(fp);
}

void VectorProcessSize(const char *path, Vector *vec)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - vector size \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }
    fscanf(fp, "%d", &(vec->n));
    fclose(fp);
}

void MatrixProcess(const char *path, Matrix *mat, int row_start, int row_end)
{
    printf("loc_row_start = %d, loc_row_end = %d\n", row_start, row_end);
    int nnz_loc = 0;

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
    fgets(buffer, MAX_SIZE, fp);

    for (int index = 0; index < mat->nnz; ++index)
    {
        int m_tmp = 0, n_tmp = 0;
        double val_tmp_re = 0., val_tmp_im = 0.;
        fscanf(fp, "%d%d%lf%lf", &m_tmp, &n_tmp, &val_tmp_re, &val_tmp_im);
        --m_tmp;
        if (m_tmp >= row_start && m_tmp < row_end)
        {
            ++nnz_loc;
        }
    }
    printf("---- local nnz = %d\n", nnz_loc);
    fclose(fp);

    // getting local row_idx, col_idx, val
    if ((mat->row_idx = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (mat->col_idx = (int *)malloc(nnz_loc * sizeof(int))) == NULL ||
        (mat->val_re = (double *)malloc(nnz_loc * sizeof(double))) == NULL ||
        (mat->val_im = (double *)malloc(nnz_loc * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'matrix data\'\n");
        exit(EXIT_FAILURE);
    }

    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix process \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    // skip first 3 lines
    fgets(buffer, MAX_SIZE, fp);
    fgets(buffer, MAX_SIZE, fp);
    fgets(buffer, MAX_SIZE, fp);

    int loc_count = 0;
    for (int index = 0; index < mat->nnz; ++index)
    {
        int m_tmp = 0, n_tmp = 0;
        double val_tmp_re = 0., val_tmp_im = 0.;
        fscanf(fp, "%d%d%lf%lf", &m_tmp, &n_tmp, &val_tmp_re, &val_tmp_im);
        --m_tmp;
        if (m_tmp >= row_start && m_tmp < row_end)
        {
            // base-1 to base-0
            mat->row_idx[loc_count] = m_tmp;
            mat->col_idx[loc_count] = n_tmp - 1;
            mat->val_re[loc_count] = val_tmp_re;
            mat->val_im[loc_count] = val_tmp_im;
            ++loc_count;
        }
    }
    fclose(fp);

    // updating nnz to local nnz
    mat->nnz = nnz_loc;
}

void VectorProcess(const char *path, Vector *vec, int row_start, int row_end)
{
    printf("loc_row_start = %d, loc_row_end = %d\n", row_start, row_end);
    double *val_tmp_re = NULL, *val_tmp_im = NULL;
    int loc_size = row_end - row_start;

    if ((val_tmp_re = (double *)malloc(vec->n * sizeof(double))) == NULL ||
        (val_tmp_im = (double *)malloc(vec->n * sizeof(double))) == NULL ||
        (vec->val_re = (double *)malloc(loc_size * sizeof(double))) == NULL ||
        (vec->val_im = (double *)malloc(loc_size * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'vector data\'\n");
        exit(EXIT_FAILURE);
    }

    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - vector data \'%s\'\n", path);
        exit(EXIT_FAILURE);
    }

    int n_tmp = 0;
    fscanf(fp, "%d", &n_tmp);
    for (int index = 0; index < n_tmp; ++index)
    {
        fscanf(fp, "%lf%lf", val_tmp_re + index, val_tmp_im + index);
    }

    fclose(fp);

    for (int index = row_start; index < row_end; ++index)
    {
        vec->val_re[index - row_start] = val_tmp_re[index];
        vec->val_im[index - row_start] = val_tmp_im[index];
    }

    // updating size to local size
    vec->n = loc_size;

    // free memory
    free(val_tmp_re);
    free(val_tmp_im);
}
