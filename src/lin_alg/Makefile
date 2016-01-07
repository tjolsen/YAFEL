include $(YAFELDIR)/common.mk

cppfiles = $(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: .depend $(ofiles)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm $(ofiles) .depend
