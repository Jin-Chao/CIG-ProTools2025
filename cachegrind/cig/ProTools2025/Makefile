SUBDIRS := Himeno/src ./PolyBenchC-4.2.1 ./Simple/src ./IRSmk/IRSmk.ori ./IRSmk/IRSmk.aos

.PHONY: all clean $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

clean:
	for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir clean; \
	done

