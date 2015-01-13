subdirs = src apps doc

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@

clean:
	for dir in $(subdirs); do $(MAKE) -C $$dir clean; done
	rm $(YAFELDIR)/lib/$(LIB)

.PHONY: all $(subdirs)