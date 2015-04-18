cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
outfile = advect2D.out


.PHONY: all clean

all: $(outfile)

include $(YAFELDIR)/common.mk

$(outfile): $(ofiles)
	$(CPP) $(ofiles) $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(outfile) *.o