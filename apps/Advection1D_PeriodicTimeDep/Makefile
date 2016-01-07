cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefile = advDiffSolve.out

.PHONY: all clean

all: $(exefile)

include $(YAFELDIR)/common.mk

$(exefile): $(ofiles)
	$(CPP) $(ofiles) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm *.o *.vtu $(exefile)
