include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)

all: .depend $(exefiles)


%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm .depend $(ofiles) $(exefiles)