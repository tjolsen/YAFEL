include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: $(ofiles)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(ofiles)