cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

include $(YAFELDIR)/common.mk

all: .depend $(ofiles)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(ofiles)

.PHONY: all clean