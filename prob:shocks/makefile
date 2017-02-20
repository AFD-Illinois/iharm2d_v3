#
#
CORE_SRC_DIR = ../core_src
#
CC      = gcc
CFLAGS  = -Wall -Ofast -fopenmp -I./
LDFLAGS = -lm -lgsl -lgslcblas
EXE     = harm

# core routines
CORE_SRC = $(wildcard $(CORE_SRC_DIR)/*.c)

# problem dependent routines
PROB_SRC = $(wildcard ./*.c)

SRC = $(CORE_SRC) $(PROB_SRC)
OBJS = $(SRC:.c=.o) 
INCS = $(wildcard ./*.h)

# makes harm
$(EXE): $(OBJS) $(INCS) makefile
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(EXE)

$(OBJS) : $(INCS) $(SRCS) makefile

.PHONY: clean
clean:
	/bin/rm -f *.o 
	/bin/rm -f $(CORE_SRC_DIR)/*.o

newrun:
	/bin/rm -rf dumps images ener.out harm.out
	/bin/rm -f *.o
	/bin/rm -f $(CORE_SRC_DIR)/*.o
	make harm

