testfiles = test_vectorTest.cpp test_matrixTest.cpp test_matmul.cpp
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
	@for f in $(test_exefiles); do echo "Running: $$f"; ./$$f; done

runbenchmarks: all
	@for f in $(benchmark_exefiles); do echo "Running: $$f"; ./$$f; done

.PHONY: all clean runtests runbenchmarks