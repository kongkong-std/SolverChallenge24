#ifndef _MMIO_HIGHLEVEL_
#define _MMIO_HIGHLEVEL_

#ifndef VALUE_TYPE
#define VALUE_TYPE double
#endif

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif
#ifndef MAT_VAL_TYPE
#define MAT_VAL_TYPE VALUE_TYPE
#endif
#include "mmio.h"
void exclusive_scan(MAT_PTR_TYPE *input, int length);

int read_mtx_header(char *filename, int *isComplex);

int mmio_allinone(int *m, int *n, int *nnz, int *isSymmetric, int *base,
                  int **csrRowPtr, int **csrColIdx, MAT_VAL_TYPE **csrVal, char *filename);

int mmio_allinone_complex(int *m, int *n, int *nnz, int *isSymmetric, int *base, int **csrRowPtr,
                          int **csrColIdx, MAT_VAL_TYPE **csrVal_re,
                          MAT_VAL_TYPE **csrVal_im, char *filename);
#endif
