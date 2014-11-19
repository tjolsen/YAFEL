include $(YAFELDIR)/common.mk
include $(YAFELDIR)/apps/appscommon.mk

output = NS_Chorin
cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)

$(output): $(ofiles)
	$(CPP) -o $@ $(ofiles) $(LFLAGS)

%.o: %.cpp
	$(CPP) $(CFLAGS) $<

clean:
	rm $(output) $(ofiles) *.vtu