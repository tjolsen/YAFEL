
cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)

include $(YAFELDIR)/common.mk


all: .depend $(exefiles)


%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) .depend