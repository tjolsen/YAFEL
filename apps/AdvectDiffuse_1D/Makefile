cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefile = advDiffSolve.out

.PHONY: all clean

all: $(exefile)

include $(YAFELDIR)/common.mk

$(exefile): $(ofiles)
	$(CXX) $(ofiles) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm *.o  $(exefile)
