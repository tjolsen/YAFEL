include $(YAFELDIR)/common.mk
include $(YAFELDIR)/apps/appscommon.mk

output = Truss2D_solver
cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

all: .depend $(output)

$(output): $(ofiles)
	$(CXX) -o $@ $(ofiles) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $<

clean:
	rm $(output) $(ofiles) *.vtu .depend
