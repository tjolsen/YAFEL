include $(YAFELDIR)/common.mk

cppfiles = $(shell ls *.cpp)
ofiles = $(cppfiles:.cpp=.o)
output = test_square_grid

all: .depend $(output)

$(output): $(ofiles)
	$(CXX) $^ $(LDFLAGS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $< 

clean:
	rm $(output) $(ofiles) *.vtu
