# makefile for check
CC = gcc
SRC = ../src

CFLAGS = -c -std=c99 -Wall -Wextra
OBJECT := $(patsubst %.c,%.o,$(wildcard $(SRC)/*.c))

CMD =
INCL =
LIB_PATH =
LIB_FLAG = -lm

LIB = $(LIB_PATH) $(LIB_FLAG)

.PHONY: all clean define_csr define_check_error

# Display progress
TOTAL_FILES := $(shell echo $(OBJECT) | wc -w)
COUNTER = 0
define increment_counter
	$(eval COUNTER=$(shell echo $$(($(COUNTER)+1))))
	echo -n "[$(COUNTER)/$(TOTAL_FILES)] Compiling $< ... "
endef

all: app_check

app_check: $(OBJECT)
	$(CC) $(CMD) $(INCL) -o app_check $(OBJECT) $(LIB)

%.o: $(SRC)/%.c
	$(increment_counter)
	$(CC) $(CFLAGS) $(CMD) $(INCL) $(DEFINES) -o "$@" "$<"
	echo "Done."

clean:
	rm -f app_check $(OBJECT)

define_csr:
	$(MAKE) all DEFINES=-DCSR_MATRIX_FILE_

define_check_error:
	$(MAKE) all DEFINES=-DCHECK_ERROR_