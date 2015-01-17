include $(YAFELDIR)/common.mk
include $(YAFELDIR)/apps/appscommon.mk

output = Truss2D_solver
cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: .depend $(output)

$(output): $(ofiles)
	$(CPP) -o $@ $(ofiles) $(LFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(output) $(ofiles) *.vtu .depend