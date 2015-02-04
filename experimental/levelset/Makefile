include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
output = test_square_grid

all: .depend $(output)

$(output): $(ofiles)
	$(CPP) $^ $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $(CFLAGS) $< 

clean:
	rm $(output) $(ofiles) *.vtu