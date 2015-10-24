#parallelize with OpenMP
useOpenMP = false;

#compile lin_alg routines in "optimized" mode:
#  disables certain bounds checks
linalg_optimized = true;

#define C++ compiler
CPP = g++

#define C compiler
CC = gcc

#define archive utility
AR = ar

#compiler optimization and linking flags
CFLAGS = -O3 -march=native -Wall
LFLAGS = -L$(YAFELDIR)/lib/ -lyafel
ARFLAGS = -ru

#output library name
LIB = libyafel.a


#=================================================================
# DON'T TOUCH ANYTHING BELOW HERE
#=================================================================

# These will always be necessary, so they're moved away from the configurable line
CFLAGS += -c -I$(YAFELDIR)/include/ -std=c++11

ifeq ($(useOpenMP), true)
	CFLAGS += -fopenmp
	LFLAGS += -fopenmp
endif

ifeq ($(linalg_optimized), true)
	CFLAGS += -D_OPTIMIZED
endif

#======================================================
# Make dependencies
#======================================================
.depend:
	$(CPP) $(shell ls *.cpp) $(CFLAGS) -MM > .depend

-include .depend
