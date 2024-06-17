#include "main.h"
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

    for (int index = 0; index < mat->nnz; ++index)
    {
        int m_tmp = 0;
        fscanf(fp, "%d%*d%*lf", &m_tmp);
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
        (mat->val = (double *)malloc(nnz_loc * sizeof(double))) == NULL)
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
        double val_tmp = 0.;
        fscanf(fp, "%d%d%lf", &m_tmp, &n_tmp, &val_tmp);
        --m_tmp;
        if (m_tmp >= row_start && m_tmp < row_end)
        {
            // base-1 to base-0
            mat->row_idx[loc_count] = m_tmp;
            mat->col_idx[loc_count] = n_tmp - 1;
            mat->val[loc_count] = val_tmp;
            ++loc_count;
        }
    }
    fclose(fp);

    // updating nnz to local nnz
    mat->nnz = nnz_loc;

#if 0
    MPI_File fh;
    MPI_Status status;
    MPI_Offset file_size;

    MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &file_size);

    char *buffer = (char *)malloc(file_size * sizeof(char));
    MPI_File_read_at_all(fh, 0, buffer, file_size, MPI_CHAR, &status);
    MPI_File_close(&fh);

    char *ptr = buffer;
    // Skip the first two lines
    ptr = strstr(ptr, "\n") + 1;
    ptr = strstr(ptr, "\n") + 1;

    int m, n, nnz;
    sscanf(ptr, "%d %d %d", &m, &n, &nnz);
    ptr = strstr(ptr, "\n") + 1;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int local_nnz = nnz / size;
    int remainder = nnz % size;
    if (rank < remainder)
    {
        local_nnz++;
    }

    mat->m = m;
    mat->n = n;
    mat->nnz = local_nnz;

    mat->row_idx = (int *)malloc(local_nnz * sizeof(int));
    mat->col_idx = (int *)malloc(local_nnz * sizeof(int));
    mat->val = (double *)malloc(local_nnz * sizeof(double));

    // Skip data for other processes and read data for this process
    for (int i = 0; i < rank * local_nnz; ++i)
    {
        ptr = strstr(ptr, "\n") + 1;
    }
    for (int index = 0; index < local_nnz; ++index)
    {
        sscanf(ptr, "%d %d %lf", &mat->row_idx[index], &mat->col_idx[index], &mat->val[index]);
        ptr = strstr(ptr, "\n") + 1;
    }

    free(buffer);
#endif
}

void VectorProcess(const char *path, Vector *vec, int row_start, int row_end)
{
    printf("loc_row_start = %d, loc_row_end = %d\n", row_start, row_end);
    double *val_tmp = NULL;
    int loc_size = row_end - row_start;

    if ((val_tmp = (double *)malloc(vec->n * sizeof(double))) == NULL ||
        (vec->val = (double *)malloc(loc_size * sizeof(double))) == NULL)
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
    for(int index = 0; index < n_tmp; ++index)
    {
        fscanf(fp, "%lf", val_tmp + index);
    }

    fclose(fp);

    for(int index = row_start; index < row_end; ++index)
    {
        vec->val[index - row_start] = val_tmp[index];
    }

    // updating size to local size
    vec->n = loc_size;

    // free memory
    free(val_tmp);

#if 0
    MPI_File fh;
    MPI_Status status;
    MPI_Offset file_size;

    MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &file_size);

    char *buffer = (char *)malloc(file_size * sizeof(char));
    MPI_File_read_at_all(fh, 0, buffer, file_size, MPI_CHAR, &status);
    MPI_File_close(&fh);

    char *ptr = buffer;

    int n;
    sscanf(ptr, "%d", &n);
    ptr = strstr(ptr, "\n") + 1;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int local_n = n / size;
    int remainder = n % size;
    if (rank < remainder)
    {
        local_n++;
    }

    vec->n = local_n;
    vec->val = (double *)malloc(local_n * sizeof(double));

    for (int i = 0; i < rank * local_n; ++i)
    {
        ptr = strstr(ptr, "\n") + 1;
    }
    for (int index = 0; index < local_n; ++index)
    {
        sscanf(ptr, "%lf", &vec->val[index]);
        ptr = strstr(ptr, "\n") + 1;
    }

    free(buffer);
#endif
}
