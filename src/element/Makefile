cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

include $(YAFELDIR)/common.mk

all: $(ofiles)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(ofiles)

.PHONY: all clean