testfiles = test_vectorTest.cpp test_matrixTest.cpp test_sparse_coo.cpp test_sparse_csr.cpp test_tensor.cpp #test_matmul.cpp 
benchmarkfiles = matmul_performance.cpp
cppfiles = $(testfiles) $(benchmarkfiles)

ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)
test_exefiles = $(testfiles:.cpp=.exe)
benchmark_exefiles = $(benchmarkfiles:.cpp=.exe)

all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) .depend

runtests: all
	@for f in $(test_exefiles); do echo "Running: $$f"; (ulimit -S -s 8192; ./$$f); done 

runbenchmarks: all
	(ulimit -a)
	@for f in $(benchmark_exefiles); do echo "Running: $$f"; (ulimit -S -s 8192; a=./$$f); done

.PHONY: all clean runtests runbenchmarks
