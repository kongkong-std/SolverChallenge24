#include "main.h"

int main(int argc, char **argv)
{
    int np = 0;
    char *dst_sol = NULL;
    for (int index = 0; index < argc; ++index)
    {
        if (strstr("-np", argv[index]))
        {
            np = atoi(argv[index + 1]);
        }
        if (strstr("-dst_sol", argv[index]))
        {
            dst_sol = argv[index + 1];
        }
    }

    Vector *answer_array = NULL;
    if ((answer_array = (Vector *)malloc(np * sizeof(Vector))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'answer array\'\n");
        exit(EXIT_FAILURE);
    }

    char **answer_file_array;
    if ((answer_file_array = (char **)malloc(np * sizeof(char *))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'answer file **\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < np; ++index)
    {
        if ((*(answer_file_array + index) = (char *)malloc(MAX_BUFFER_SIZE * sizeof(char))) == NULL)
        {
            fprintf(stderr, "Memory allocation failed - \'answer file *\'\n");
            exit(EXIT_FAILURE);
        }
    }

    for (int index = 0; index < np; ++index)
    {
        sprintf(answer_file_array[index], "./answer_file/answer_x_%d.txt", index);
        AnswerXFileProcess(answer_file_array[index], &(answer_array[index]));
    }

    int sol_n = 0;
    for (int index = 0; index < np; ++index)
    {
        sol_n += answer_array[index].n;
    }
    printf(">>>> dimension of solution: %d\n", sol_n);

    Vector sol_vec;
    sol_vec.n = sol_n;
#ifdef COMPLEX_VECTOR_
    if ((sol_vec.row_idx = (int *)malloc(sol_n * sizeof(int))) == NULL ||
        (sol_vec.val_re = (double *)malloc(sol_n * sizeof(double))) == NULL ||
        (sol_vec.val_im = (double *)malloc(sol_n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'solution vector file\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < np; ++index)
    {
        for (int index_j = 0; index_j < answer_array[index].n; ++index_j)
        {
            sol_vec.row_idx[answer_array[index].row_idx[index_j]] = answer_array[index].row_idx[index_j];
            sol_vec.val_re[answer_array[index].row_idx[index_j]] = answer_array[index].val_re[index_j];
            sol_vec.val_im[answer_array[index].row_idx[index_j]] = answer_array[index].val_im[index_j];
        }
    }
#elif defined REAL_VECTOR_
    if ((sol_vec.row_idx = (int *)malloc(sol_n * sizeof(int))) == NULL ||
        (sol_vec.val = (double *)malloc(sol_n * sizeof(double))) == NULL)
    {
        fprintf(stderr, "Memory allocation failed - \'solution vector file\'\n");
        exit(EXIT_FAILURE);
    }

    for (int index = 0; index < np; ++index)
    {
        for (int index_j = 0; index_j < answer_array[index].n; ++index_j)
        {
            sol_vec.row_idx[answer_array[index].row_idx[index_j]] = answer_array[index].row_idx[index_j];
            sol_vec.val[answer_array[index].row_idx[index_j]] = answer_array[index].val[index_j];
        }
    }
#endif

    FILE *fp = NULL;

    if ((fp = fopen(dst_sol, "wb")) == NULL)
    {
        fprintf(stderr, "Cannot open file - \'%s\'\n", dst_sol);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "%d\n", sol_vec.n);
#ifdef COMPLEX_VECTOR_
    for (int index = 0; index < sol_vec.n; ++index)
    {
        // fprintf(fp, "%d\t%021.16le\t%021.16le\n", sol_vec.row_idx[index],
        fprintf(fp, "%021.16le\t%021.16le\n", sol_vec.val_re[index], sol_vec.val_im[index]);
    }
#elif defined REAL_VECTOR_
    for (int index = 0; index < sol_vec.n; ++index)
    {
        // fprintf(fp, "%d\t%021.16le\t%021.16le\n", sol_vec.row_idx[index],
        fprintf(fp, "%021.16le\n", sol_vec.val[index]);
    }
#endif
    fclose(fp);

    // free memory
    free(sol_vec.row_idx);
#ifdef COMPLEX_VECTOR_
    free(sol_vec.val_re);
    free(sol_vec.val_im);
#elif defined REAL_VECTOR_
    free(sol_vec.val);
#endif
    for (int index = 0; index < np; ++index)
    {
        free(*(answer_file_array + index));
    }
    free(answer_file_array);
    for (int index = 0; index < np; ++index)
    {
        free(answer_array[index].row_idx);
#ifdef COMPLEX_VECTOR_
        free(answer_array[index].val_re);
        free(answer_array[index].val_im);
#elif defined REAL_VECTOR_
        free(answer_array[index].val);
#endif
    }
    free(answer_array);

    return 0;
}

// command
/*
 * ./a.out -np <int> -dst_sol <path/to/solution>
 */