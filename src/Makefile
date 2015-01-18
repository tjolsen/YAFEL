SUBDIRS = lin_alg mesh element output utils
SUBDIRTARGETS= $(SUBDIRS:%=%.tmp) # workaround to get header-file dependencies to work correctly

include $(YAFELDIR)/common.mk

all: $(LIB) stage
	echo $(SUBDIRTARGETS)

$(LIB): $(SUBDIRTARGETS)
	$(AR) $(ARFLAGS) $@ $(shell find . | grep .o$$)

stage: $(LIB)
	cp $(LIB) $(YAFELDIR)/lib/

$(SUBDIRTARGETS):
	$(MAKE) -C $(@:.tmp=) all

clean:
	for dir in $(SUBDIRS); do $(MAKE) -C $$dir clean; done 
	rm $(LIB)

.PHONY: all stage clean $(SUBDIRTARGETS)