# make file for prodecomp2

#gcc -g -pg -Wall -pedantic -std=c99 -Ofast readFT2.c decompose2.c fastnnls.c midmatrix3.c setRandom.c startShape.c getParam.c rmsd.c util.c readInput.c readPeaklist.c util2.c correlation.c  -fopenmp main2.c -L/usr/lib/openblas-base -lgsl -lopenblas -lm -o main 

# gcc -g -pg -O3 -Wall -pedantic -std=c99 readFT2.c decompose2.c fastnnls.c midmatrix3.c setRandom.c startShape.c getParam.c rmsd.c util.c readInput.c readPeaklist.c util2.c correlation.c -fopenmp main2.c -lgsl -lcblas -latlas -lm -o main

CC=gcc
CFLAGS=-g -pg -Wall -pedantic -std=c99 -Ofast -fopenmp
INCLUDES=common.h
LFLAGS=-L/usr/lib/openblas-base -lgsl -lopenblas -lm
#LFLAGS=-lgsl -lcblas -latlas -lm
SRCS=readFT2.c decompose3.c fastnnls3.c midmatrix2.c setRandom.c startShape.c getParam.c rmsd.c util.c readInput.c readPeaklist.c util2.c correlation.c main3.c
OBJS=$(SRCS:.c=.o)
TARGET=main

all:    $(TARGET)
	@echo  compiling main

$(TARGET): $(OBJS) 
	$(CC) $(CFLAGS) $(SRCS) $(LFLAGS) -o $(TARGET)

clean:	
	$(RM) *.o
	@echo  removing object files	
