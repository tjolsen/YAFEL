cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = surface_area.exe
ofiles := $(filter-out $(exefiles:.exe=.o), $(ofiles))



all: .depend $(exefiles)

include $(YAFELDIR)/common.mk
CXXFLAGS += -D_PARALLEL_ILU

%.exe: %.o $(ofiles)
	$(CXX) $< $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm $(ofiles) $(exefiles) *.vtu *.pvd $(exefiles:.exe=.o)

run: all
	./$(exefile)

.PHONY: all clean run
