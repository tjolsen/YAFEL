SUBDIRS = lin_alg mesh element output

include $(YAFELDIR)/common.mk

all: $(LIB) stage

$(LIB): $(SUBDIRS)
	$(AR) $(ARFLAGS) $@ $(shell find . | grep .o$$)
	ranlib $@

stage: $(LIB)
	cp $(LIB) $(YAFELDIR)/lib/

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done 
	rm $(LIB)

.PHONY: all stage clean $(SUBDIRS)