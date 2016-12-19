#parallelize with OpenMP
useOpenMP = false

#compile lin_alg routines in "optimized" mode:
#  disables certain bounds checks
linalg_optimized = true

# Use Parallel Matrix multiplication algorithm
parallel_matmul = true

#define C++ compiler
ifndef CXX
 export CXX = g++
endif

#define C compiler
ifndef CC
 export CC = gcc
endif

#define archive utility
AR = ar

#compiler optimization and linking flags
CXXFLAGS = -O3 -Wall -Werror -Wextra #-fmax-errors=5 
LDFLAGS = -L$(YAFELDIR)/lib/ -lyafel
ARFLAGS = -ru

#output library name
LIB = libyafel.a


#=================================================================
# DON'T TOUCH ANYTHING BELOW HERE
#=================================================================

# These will always be necessary, so they're moved away from the configurable line
CXXFLAGS += -c -I. -I$(YAFELDIR)/include -std=c++11 -march=native -mtune=native -funroll-loops -mavx -ffast-math

ifeq ($(useOpenMP), true)
	CXXFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

ifeq ($(linalg_optimized), true)
	CXXFLAGS += -D_OPTIMIZED
endif

ifeq ($(parallel_matmul), true)
	CXXFLAGS += -pthread -D_YAFEL_PARALLEL_MATMUL
	LDFLAGS += -pthread
endif
#======================================================
# Make dependencies
#======================================================
.depend:
	$(CXX) $(wildcard *.cpp) $(CXXFLAGS) -MM > .depend

-include .depend
