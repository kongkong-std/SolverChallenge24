#ifndef FILE_PROCESS_H_
#define FILE_PROCESS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// function prototype
void RealCOO2CSRMatrixFileProcess(const char * /*path to file*/, int * /*row*/, int * /*column*/, int * /*nnz*/,
                                  int ** /*row pointer*/, int ** /*column index pointer*/, double ** /*value*/);

void ComplexCOO2CSRMatrixFileProcess(const char * /*path to file*/, int * /*row*/, int * /*column*/, int * /*nnz*/,
                                     int ** /*row pointer*/, int ** /*column index pointer*/,
                                     double ** /*re value*/, double ** /*im value*/);

#endif