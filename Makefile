CC = gcc
CFLAGS = -lm -Ofast -I$(DIR_INC)
DFLAGS = -Wall -Wextra -W -g -O0 -lm -I$(DIR_INC)

DIR_SRC = src
DIR_BIN = bin
DIR_INC = include

_DEPS = functions.h
DEPS = $(patsubst %,$(DIR_INC)/%,$(_DEPS))

_OBJ = main.o functions.o
OBJ = $(patsubst %,$(DIR_BIN)/%,$(_OBJ))

all : check_dir $(DIR_BIN)/main

$(DIR_BIN)/main : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

$(DIR_BIN)/%.o : $(DIR_SRC)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY : clean check_dir

clean :
	rm -f $(DIR_BIN)/*.o $(DIR_BIN)/main

check_dir :
	if [ ! -d $(DIR_BIN) ]; then mkdir $(DIR_BIN); fi
