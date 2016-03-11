cppfiles = main.cpp LinearAdvection.cpp
ofiles = $(cppfiles:.cpp=.o)
exefile = linear_advection.exe

all: .depend $(exefile)

include $(YAFELDIR)/common.mk

$(exefile): $(ofiles)
	$(CXX) $(ofiles) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm $(ofiles) $(exefile) *.vtu *.pvd *.csv

run: all
	./$(exefile)

.PHONY: all clean run
