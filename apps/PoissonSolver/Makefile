include $(YAFELDIR)/common.mk
#include $(YAFELDIR)/apps/appscommon.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
output = poissonSolver

all: .depend $(output)

$(output): $(ofiles)
	$(CXX) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< 

clean:
	rm $(output) $(ofiles) *.vtu
