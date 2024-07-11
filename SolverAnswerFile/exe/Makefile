# makefile for check
CC = gcc
SRC = ../src

CFLAGS = -c -std=c99
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

CMD =
INCL =
LIB_PATH =
LIB_FLAG =

LIB = $(LIB_PATH) $(LIB_FLAG)

DEFINES =

.PHONY: all clean define_complex_vector define_real_vector

all: app_answer_file

app_answer_file: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_answer_file $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) $(DEFINES) "$<" -o "$@"

clean:
	rm -f app_answer_file $(OBJECT)

define_complex_vector:
	$(MAKE) all DEFINES=-DCOMPLEX_VECTOR_

define_real_vector:
	$(MAKE) all DEFINES=-DREAL_VECTOR_