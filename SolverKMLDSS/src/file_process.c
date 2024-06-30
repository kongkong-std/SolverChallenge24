#include "file_process.h"
#define MAX_BUFFER_SIZE 4096

void RealRHSFileProcess(const char *path, double **rhs)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - rhs file\n");
        exit(EXIT_FAILURE);
    }

    char buffer[MAX_BUFFER_SIZE];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    int n = 0;
    fscanf(fp, "%d", &n);

    if ((*rhs = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - rhs data\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        fscanf(fp, "%lf", *rhs + index);
    }

    fclose(fp);
}

void ComplexRHSFileProcess(const char *path, double **re_rhs, double **im_rhs)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - rhs file\n");
        exit(EXIT_FAILURE);
    }

    char buffer[MAX_BUFFER_SIZE];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    int n = 0;
    fscanf(fp, "%d", &n);

    if ((*re_rhs = (double *)malloc(n * sizeof(double))) == NULL ||
        (*im_rhs = (double *)malloc(n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - rhs data\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < n; ++index)
    {
        fscanf(fp, "%lf%lf", *re_rhs + index, *im_rhs + index);
    }

    fclose(fp);
}

void RealCOO2CSRMatrixFileProcess(const char *path, int *m, int *n, int *nnz,
                                  int **csr_row_ptr, int **csr_col_idx, double **csr_val)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix file\n");
        exit(EXIT_FAILURE);
    }

    char buffer[MAX_BUFFER_SIZE];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    int m_mat = 0, n_mat = 0, nnz_mat = 0;
    fscanf(fp, "%d%d%d", &m_mat, &n_mat, &nnz_mat);
    *m = m_mat;
    *n = n_mat;
    *nnz = nnz_mat;

    int *row_idx = NULL, *col_idx = NULL;
    double *val = NULL;
    if ((row_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (col_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (val = (double *)malloc(nnz_mat * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - coo matrix data\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < nnz_mat; ++index)
    {
        fscanf(fp, "%d%d%lf", row_idx + index, col_idx + index, val + index);
        row_idx[index] = row_idx[index] - 1;
        col_idx[index] = col_idx[index] - 1;
    }

    fclose(fp);

    if ((*csr_row_ptr = (int *)malloc((n_mat + 1) * sizeof(int))) == NULL ||
        (*csr_col_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (*csr_val = (double *)malloc(nnz_mat * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocatio failed - csr matrix data\n");
        exit(EXIT_FAILURE);
    }

    memset(*csr_row_ptr, 0, (n_mat + 1) * sizeof(int));

    // Count the number of entries in each row
    for (int index = 0; index < nnz_mat; ++index)
    {
        (*csr_row_ptr)[row_idx[index] + 1] = (*csr_row_ptr)[row_idx[index] + 1] + 1;
    }

    // Cumulative sum to get row_ptr
    for (int index = 0; index < m_mat; ++index)
    {
        (*csr_row_ptr)[index + 1] = (*csr_row_ptr)[index + 1] + (*csr_row_ptr)[index];
    }

    // Fill csr_col_idx and csr_val
    int *temp_row_ptr = (int *)malloc((m_mat + 1) * sizeof(int));
    memcpy(temp_row_ptr, *csr_row_ptr, (m_mat + 1) * sizeof(int));

    for (int index = 0; index < nnz_mat; ++index)
    {
        int row = row_idx[index];
        int dest = temp_row_ptr[row];

        (*csr_col_idx)[dest] = col_idx[index];
        (*csr_val)[dest] = val[index];

        temp_row_ptr[row]++;
    }

    free(temp_row_ptr);

    // free memory
    free(row_idx);
    free(col_idx);
    free(val);
}

void ComplexCOO2CSRMatrixFileProcess(const char *path, int *m, int *n, int *nnz,
                                     int **csr_row_ptr, int **csr_col_idx,
                                     double **csr_val_re, double **csr_val_im)
{
    FILE *fp = NULL;
    if ((fp = fopen(path, "rb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - matrix file\n");
        exit(EXIT_FAILURE);
    }

    char buffer[MAX_BUFFER_SIZE];
    while (fgets(buffer, sizeof(buffer), fp) != NULL)
    {
        if (buffer[0] != '%')
        {
            // Move file pointer back to the beginning of this line
            fseek(fp, -strlen(buffer), SEEK_CUR);
            break;
        }
    }

    int m_mat = 0, n_mat = 0, nnz_mat = 0;
    fscanf(fp, "%d%d%d", &m_mat, &n_mat, &nnz_mat);
    *m = m_mat;
    *n = n_mat;
    *nnz = nnz_mat;

    int *row_idx = NULL, *col_idx = NULL;
    double *val_re = NULL, *val_im = NULL;
    if ((row_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (col_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (val_re = (double *)malloc(nnz_mat * sizeof(double))) == NULL ||
        (val_im = (double *)malloc(nnz_mat * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - coo matrix data\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < nnz_mat; ++index)
    {
        fscanf(fp, "%d%d%lf%lf", row_idx + index, col_idx + index, val_re + index, val_im + index);
        row_idx[index] = row_idx[index] - 1;
        col_idx[index] = col_idx[index] - 1;
    }

    fclose(fp);

    if ((*csr_row_ptr = (int *)malloc((n_mat + 1) * sizeof(int))) == NULL ||
        (*csr_col_idx = (int *)malloc(nnz_mat * sizeof(int))) == NULL ||
        (*csr_val_re = (double *)malloc(nnz_mat * sizeof(double))) == NULL ||
        (*csr_val_im = (double *)malloc(nnz_mat * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocatio failed - csr matrix data\n");
        exit(EXIT_FAILURE);
    }

    memset(*csr_row_ptr, 0, (n_mat + 1) * sizeof(int));

    // Count the number of entries in each row
    for (int index = 0; index < nnz_mat; ++index)
    {
        (*csr_row_ptr)[row_idx[index] + 1] = (*csr_row_ptr)[row_idx[index] + 1] + 1;
    }

    // Cumulative sum to get row_ptr
    for (int index = 0; index < m_mat; ++index)
    {
        (*csr_row_ptr)[index + 1] = (*csr_row_ptr)[index + 1] + (*csr_row_ptr)[index];
    }

    // Fill csr_col_idx and csr_val
    int *temp_row_ptr = (int *)malloc((m_mat + 1) * sizeof(int));
    memcpy(temp_row_ptr, *csr_row_ptr, (m_mat + 1) * sizeof(int));

    for (int index = 0; index < nnz_mat; ++index)
    {
        int row = row_idx[index];
        int dest = temp_row_ptr[row];

        (*csr_col_idx)[dest] = col_idx[index];
        (*csr_val_re)[dest] = val_re[index];
        (*csr_val_im)[dest] = val_im[index];

        temp_row_ptr[row]++;
    }

    free(temp_row_ptr);

    // free memory
    free(row_idx);
    free(col_idx);
    free(val_re);
    free(val_im);
}