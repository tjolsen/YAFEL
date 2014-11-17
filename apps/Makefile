subdirs = NS_1

include $(YAFELDIR)/common.mk

all: $(subdirs)

$(subdirs):
	$(MAKE) -C $@

clean:
	for dir in $(subdirs); do $(MAKE) -C $$dir clean; done

.PHONY: all clean $(subdirs)