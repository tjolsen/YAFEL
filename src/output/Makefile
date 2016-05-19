cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)

include $(YAFELDIR)/common.mk

all: .depend $(ofiles)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm $(ofiles) .depend

.PHONY: all clean
