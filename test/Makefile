cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)

all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) .depend