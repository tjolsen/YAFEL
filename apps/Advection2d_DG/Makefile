cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefile = linear_advection.exe

all: .depend $(exefile)

include $(YAFELDIR)/common.mk

$(exefile): $(ofiles)
	$(CXX) $< $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm $(ofiles) $(exefile) *.vtu *.pvd

run: all
	./$(exefile)

.PHONY: all clean run