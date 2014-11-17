subdirs = src #test apps

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	cd $@ && $(MAKE)

clean:
	for dir in $(subdirs); do cd $(YAFELDIR)/$$dir && $(MAKE) clean; done
	rm $(YAFELDIR)/lib/$(LIB)

.PHONY: all $(subdirs)