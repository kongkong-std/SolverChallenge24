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

.PHONY: all clean

all: app_answer_file

app_answer_file: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_answer_file $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm app_answer_file $(OBJECT)