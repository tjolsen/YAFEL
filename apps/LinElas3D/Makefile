cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefile = linElas3D.out

.PHONY: all clean

all: $(exefile)

include $(YAFELDIR)/common.mk

$(exefile): $(ofiles)
	$(CPP) $(ofiles) $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm *.o *.vtu $(exefile)
