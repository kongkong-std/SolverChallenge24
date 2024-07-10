# makefile for kml solver
# module purge
# module use /workspace/public/software/modules
# module load compilers/gcc/kunpenggcc/10.3.1/gcc10.3.1
# module load kml-dss/kmldss/kmldss-gcc
# module load kml-dss/openblas/openblas-gcc
# export OMP_PROC_BIND=close
CC = gcc
SRC = ../src

CFLAGS = -c -std=c99 -fopenmp -DKML_DSS_IR_
#CFLAGS = -c -std=c99 -fopenmp
CMD = -Wall -g -O3

INCL =

LIB_PATH =

LIB_FLAG = -lksolver -lkservice -lopenblas -lm

LIB = $(LIB_PATH) $(LIB_FLAG)
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

.PHONY: all clean

all: app_kml_exe

app_kml_exe: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_kml_exe $(OBJECT) $(LIB) -fopenmp

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm app_kml_exe $(OBJECT)