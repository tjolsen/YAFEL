subdirs = src apps test doc

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ all

clean:
	for dir in $(subdirs); do $(MAKE) -C $$dir clean; done
	rm $(YAFELDIR)/lib/$(LIB)

runtests: test
	make -C test runtests

.PHONY: all $(subdirs) runtests
