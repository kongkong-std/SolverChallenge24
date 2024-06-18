# makefile for complex petsc binary file
CC = /home/export/base/ycsc_zhaol/qinzk/online1/software/mpich-3.2.1/bin/mpicc
SRC = ../src

CFLAGS = -c -std=c99
CMD = -fPIC -Wall -Wwrite-strings -Wno-unknown-pragmas \
-Wno-lto-type-mismatch -fstack-protector -fvisibility=hidden -g3 -O3

INCL = -I/home/export/base/ycsc_zhaol/qinzk/online1/software/petsc-3.19.2/complex/include \
-I/home/export/base/ycsc_zhaol/qinzk/online1/software/mpich-3.2.1/include

LIB_PATH = -L/home/export/base/ycsc_zhaol/qinzk/online1/software/petsc-3.19.2/complex/lib \
-Wl,-rpath,/home/export/base/ycsc_zhaol/qinzk/online1/software/lapack-3.11.0 \
-L/home/export/base/ycsc_zhaol/qinzk/online1/software/lapack-3.11.0 \
-Wl,-rpath,/home/export/base/ycsc_zhaol/qinzk/online1/software/BLAS-3.11.0 \
-L/home/export/base/ycsc_zhaol/qinzk/online1/software/BLAS-3.11.0 \
-L/home/export/base/ycsc_zhaol/qinzk/online1/software/mpich-3.2.1/lib \
-L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 \
-Wl,-rpath,/home/export/base/ycsc_zhaol/qinzk/online1/software/mpich-3.2.1/lib

LIB_FLAG = -lpetsc -llapack -lblas -lpthread -lm -lX11 -ldl -lmpifort -lmpi -lgfortran \
-lm -lgfortran -lm -lgcc_s -lquadmath -ldl

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_petsc_bin

app_petsc_bin: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_petsc_bin $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm app_petsc_bin $(OBJECT)