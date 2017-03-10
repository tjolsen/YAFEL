cppfile = new_tensors.cpp
ofile = $(cppfile:.cpp=.o)
exefile = $(cppfile:.cpp=.exe)

CXX=clang++

$(exefile): $(ofile)
	$(CXX) $< -o $@

include $(YAFELDIR)/common.mk


%.o: %.cpp
	$(CXX) $< $(CXXFLAGS)

clean:
	rm $(exefile) $(ofile) .depend
