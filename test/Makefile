testfiles = test_vectorTest.cpp test_matrixTest.cpp test_sparse_coo.cpp test_sparse_csr.cpp 
testfiles += test_sparse_bcsr.cpp test_tensor.cpp test_csr_spmv.cpp test_bcsr_spmv.cpp 
testfiles += test_matmul.cpp test_DualNumber.cpp test_cg_solve.cpp test_bicgstab_solve.cpp test_ILU.cpp

benchmarkfiles = matmul_performance.cpp

cppfiles = $(testfiles) $(benchmarkfiles)

ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)
test_exefiles = $(testfiles:.cpp=.exe)
benchmark_exefiles = $(benchmarkfiles:.cpp=.exe)

all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

COVERAGE_CFLAGS += $(CFLAGS) -O0 -fprofile-arcs -ftest-coverage
COVERAGE_CXXFLAGS += $(CXXFLAGS) -O0 -fprofile-arcs -ftest-coverage
COVERAGE_LDFLAGS += $(LDFLAGS) -lgcov

test_%.exe: test_%.o
	$(CXX) $< $(COVERAGE_LDFLAGS) -o $@

test_%.o: test_%.cpp
	$(CXX) $< $(COVERAGE_CXXFLAGS)

%.exe: %.o
	$(CXX) $< $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) *.gcda *.gcno .depend coverage.info

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
