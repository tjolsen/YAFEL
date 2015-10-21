cppfiles = test_vectorTest.cpp test_matrixTest.cpp #$(wildcard *.cpp)
ofiles = $(cppfiles:.cpp=.o)
exefiles = $(cppfiles:.cpp=.exe)

all: .depend $(exefiles)

include $(YAFELDIR)/common.mk

%.exe: %.o
	$(CPP) $< $(LFLAGS) -o $@

%.o: %.cpp
	$(CPP) $< $(CFLAGS)

clean:
	rm *.o *.vtu *.pvd $(exefiles) .depend

runtests: all
	@for f in $(exefiles); do echo "Running: $$f"; ./$$f; done


.PHONY: all clean runtests