cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = laplace_solver.exe laplace_solver_instrumented.exe
ofiles := $(filter-out $(exefiles:.exe=.o), $(ofiles))


all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

%.exe: %.o $(ofiles)
	$(CXX) $< $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm $(ofiles) $(exefile) *.vtu *.pvd

run: all
	./$(exefile)

.PHONY: all clean run
