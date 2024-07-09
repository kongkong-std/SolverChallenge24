# makefile for check
CC = gcc
SRC = ../src

CFLAGS = -c
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

CMD =
INCL =
LIB_PATH =
LIB_FLAG = -lm

LIB = $(LIB_PATH) $(LIB_FLAG)

.PHONY: all clean

all: app_check

app_check: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_check $(OBJECT) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) $(CMD) $(INCL) "$<" -o "$@"

clean:
	rm app_check $(OBJECT)