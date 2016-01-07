cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
outfile = advect2D.out


.PHONY: all clean

all: $(outfile)

include $(YAFELDIR)/common.mk

$(outfile): $(ofiles)
	$(CXX) $(ofiles) $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm $(outfile) *.o
