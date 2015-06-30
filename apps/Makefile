subdirs = Truss2D_static PoissonSolver AdvectDiffuse_1D Advection1D_PeriodicTimeDep

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ all

clean:
	for dir in $(subdirs); do $(MAKE) -C $$dir clean; done

.PHONY: all clean $(subdirs)