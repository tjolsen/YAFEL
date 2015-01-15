subdirs = NS_1 Truss2D_static PoissonSolver

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@

clean:
	for dir in $(subdirs); do $(MAKE) -C $$dir clean; done

.PHONY: all clean $(subdirs)