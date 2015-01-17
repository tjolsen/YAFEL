include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: .depend $(ofiles)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(ofiles) .depend