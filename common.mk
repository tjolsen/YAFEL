#parallelize with OpenMP
useOpenMP = false

#compile lin_alg routines in "optimized" mode:
#  disables certain bounds checks
linalg_optimized = true

# Use Parallel Matrix multiplication algorithm
parallel_matmul = true

#define C++ compiler
CXX = g++

#define C compiler
CC = gcc

#define archive utility
AR = ar

#compiler optimization and linking flags
CXXFLAGS = -O3 -Wall -Werror -Wextra
LDFLAGS = #-L$(YAFELDIR)/lib/ -lyafel 
ARFLAGS = -ru

#output library name
LIB = libyafel.a


#=================================================================
# DON'T TOUCH ANYTHING BELOW HERE
#=================================================================

# These will always be necessary, so they're moved away from the configurable line
CXXFLAGS += -c -I$(YAFELDIR)/include -std=c++11 -march=native -mtune=native -funroll-loops

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
	$(CXX) $(shell ls *.cpp) $(CXXFLAGS) -MM > .depend

-include .depend
