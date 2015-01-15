include $(YAFELDIR)/common.mk
#include $(YAFELDIR)/apps/appscommon.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
output = poissonSolver

$(output): $(ofiles)
	$(CPP) $^ $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $(CFLAGS) $< 

clean:
	rm $(output) $(ofiles) *.vtu