include $(HOME)/Map/Makefile.common

ifneq (,$(shell exists.pl ./Makefile.common))
   include ./Makefile.common
endif

# You can overwrite the following variables in Makefile.common located
# in the directory make is executed from.

# CHILD_DIRS lists the subdirectories on which to recurse.
ifndef CHILD_DIRS
   CHILD_DIRS = $(shell /bin/ls -p | grep "/" | sed 's/\///' | grep -v Compendium | transpose.pl -q)
endif

CHILD_DIR1 = $(word 1, $(CHILD_DIRS))

# MAKE_TEMPLATES lists the template files that each child subdirectory
# will use as its own Makefile.
ifndef MAKE_TEMPLATES
   ifneq (,$(shell exists.pl $(MAP_TEMPLATES)/Make/$(CHILD_DIR1).mak))
      MAKE_TEMPLATES = $(foreach d, $(CHILD_DIRS), $(MAP_TEMPLATES)/Make/(d).mak)
   else
      ifneq (,$(shell exists.pl $(wildcard Lib/*.mak)))
         MAKE_TEMPLATES = $(foreach d, $(CHILD_DIRS), $(word 1, $(wildcard $(PWD)/Lib/*.mak)))
      endif
   endif
endif

indices = $(shell echo $(CHILD_DIRS) | transpose.pl -q -d ' ' | lin.pl | cut -f 1)

all:
	for child in $(CHILD_DIRS); do \
	   cd $$child && $(MAKE) ;\
	   cd .. ;\
	done

make:
	$(MKDIR) $(CHILD_DIRS); \
	$(foreach i, $(indices), \
	   cd $(word $(i), $(CHILD_DIRS)); \
	   $(RM) Makefile; \
	   $(LINK) $(word $(i), $(MAKE_TEMPLATES)) Makefile; \
	   $(MAKE) make; \
	   cd ..; ) \

make_clean:
	$(RMDIR) $(CHILD_DIRS); \

clean remote remote_clean eval_dir eval_dir_clean eval eval_clean:
	for child in $(CHILD_DIRS); do \
	   cd $$child && $(MAKE) $@ ;\
	   cd .. ;\
	done
# Note: the rules above all do the same thing, but $@ is different each time.
# $@ is always "the name of this rule." So when we run "make remote," $@ is the string "remote"

eval_clean:
	$(foreach child, $(CHILD_DIRS), \
	   cd $(child); \
	   $(MAKE) eval_clean; \
	   cd ..; )
maps: $(RUN_META)
	bind.pl $(RUN_META) | create_maps.pl -

maps_clean:
	$(RMDIR) $(CHILD_DIRS)

child_runs:
	$(foreach child, $(CHILD_DIRS), \
	   cp $(MAP_TEMPLATES)/Runs/$(RUN_META) $(child); ) \

compendium:
	cd Compendium; $(MAKE)

compendium_clean:
	cd Compendium; $(MAKE) clean


