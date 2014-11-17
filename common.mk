#parallelize with OpenMP
useOpenMP = true

#compile lin_alg routines in "optimized" mode:
#  disables certain bounds checks
linalg_optimized = true

#define C++ compiler
CPP = g++

#define C compiler
CC = gcc

#define archive utility
AR = ar

#compiling and linking flags
CFLAGS = -O3 -c -march=native -Wall -I$(YAFELDIR)/include/
LFLAGS = -L$(YAFELDIR)/lib/ -lyafel
ARFLAGS = -ru

#output library name
LIB = libyafel.a

ifeq ($(useOpenMP), true)
	CFLAGS += -fopenmp
	LFLAGS += -fopenmp
endif

ifeq ($(linalg_optimized), true)
	CFLAGS += -D_OPTIMIZED
endif
