cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

include $(YAFELDIR)/common.mk

all: .depend $(ofiles)

%.o: %.cpp
	$(CC) $(CFLAGS) $<

clean:
	rm $(ofiles) .depend

.PHONY: all clean