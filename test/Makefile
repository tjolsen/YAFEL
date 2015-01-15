include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)

all: $(exefiles)


%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles)