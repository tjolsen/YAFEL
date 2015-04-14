include $(YAFELDIR)/common.mk

cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: .depend $(ofiles)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(ofiles) .depend