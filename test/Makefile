testfiles = test_vectorTest.cpp test_matrixTest.cpp test_sparse_coo.cpp test_sparse_csr.cpp test_tensor.cpp test_csr_spmv.cpp test_matmul.cpp 
benchmarkfiles = matmul_performance.cpp
cppfiles = $(testfiles) $(benchmarkfiles)

ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)
test_exefiles = $(testfiles:.cpp=.exe)
benchmark_exefiles = $(benchmarkfiles:.cpp=.exe)

all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

COVERAGE_CFLAGS = $(CFLAGS) -O0 -fprofile-arcs -ftest-coverage
COVERAGE_CXXFLAGS = $(CFLAGS) -O0 -fprofile-arcs -ftest-coverage
COVERAGE_LDFLAGS = $(LFLAGS) -lgcov

test_%.exe: test_%.o
	$(CPP) $< $(COVERAGE_LDFLAGS) -o $@

test_%.o: test_%.cpp
	$(CPP) $< $(COVERAGE_CXXFLAGS)

%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) .depend

runtests: all
	@for f in $(test_exefiles); do echo "Running: $$f"; (ulimit -S -s 8192; ./$$f); done

coverage: all
	@for f in $(test_exefiles); do echo "Running: $$f"; (ulimit -S -s 8192; ./$$f); done
	lcov --capture --directory . --output-file coverage.info
	lcov --remove coverage.info "/usr/*" -o coverage.info
	genhtml coverage.info --output-directory outhtml

runbenchmarks: all
	(ulimit -a)
	@for f in $(benchmark_exefiles); do echo "Running: $$f"; (ulimit -S -s 8192; a=./$$f); done

.PHONY: all clean runtests runbenchmarks
