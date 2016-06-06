default: all

SHELL:=/bin/bash

include $(MAPDIR)/Templates/Make/common.mak
# If a variable called CHILDREN is defined, including this
# makefile in another makefile will provide an automatic
# way of setting up child subdirectories that all link
# to the same Makefile.  The default child makefile
# is called ../Lib/child.mak.  This can be reset by defining
# a variable called CHILD_MAKEFILE.

# Below:  ?=  sets these variables to the "default value" if they were
# not ALREADY set. If they have a value already, then ?= is nothing.
# By default, these are all blank.

# CHILDREN is a list of subdirectories that we create, do stuff in, then make the targets.
# Note that we do things in the CHILDREN directories *before* we return and make the targets
CHILDREN ?=

# "early_targets" are the files we generate BEFORE running the children
early_targets ?=

# Targets are files that we could generate with a "make (this target)" command. Note that these are generated AFTER we handle the CHILDREN
targets ?=
# side_effects are files that get generated when we're making the targets. We keep track of them here so we'll know to delete them when a user runs "make clean". But we never "make (side_effect)" any of them.
side_effects ?=
URL ?=

files_to_clean = $(early_targets) $(targets) $(side_effects)
things_to_back_up = $(early_targets) $(targets)


# Below: Set the default CHILD_MAKEFILE in case none is defined
CHILD_MAKEFILE ?= ../Lib/child.mak
ifeq (,$(CHILD_MAKEFILE))
   # Set the default CHILD_MAKEFILE in case none is defined in the file that includes quick.mak...
   CHILD_MAKEFILE = ../Lib/child.mak
endif

# === DONE setting default blank values ===

# "PHONY" targets are all rules that aren't also the name of a file.
# e.g. the "all" rule does not generate a file named "all," so it's a phony target
# (if you don't specify these, then a file named "all" can sabotage "make all" from running.
# You will get an error like "make: "all" is already up-to-date")
.PHONY: all clean back children subdirs make remote remote_clean back_clean make_clean what_settings

url_ndx        = $(shell seq 1 $(words $(URL)))
SUBDIRS        = $(CHILDREN)

# The subdir where non-local files are stored.
REMOTE_DIR ?= Remote

# The URL-related information is for automatically downloading "Remote" data to our local
# drive, then using it in some makefile step.
ifneq (,$(URL))
   REMOTE_TARGET  ?= $(foreach u, $(URL), $(addprefix Remote/, $(notdir $(u))))
   REMOTE_TARGET1 ?= $(REMOTE_DIR)/$(notdir $(word 1, $(URL)))
else
   REMOTE_TARGET  ?= $(wildcard $(REMOTE_DIR)/*)
endif

# To be effectual, need to define data.tab (i.e. fill
# in the data.tab rule at the bottom).

ifneq (,$(CHILDREN))  # NESTED IF 1 (check for non-blank CHILDREN)

delete:
	rm -rf $(CHILDREN) $(files_to_clean)

make_clean: # dangerous! kills ALL the children directories
	rm -rf $(CHILDREN);


ifneq (,$(targets))  # NESTED IF 2 (check for non-blank targets)
# If the targets aren't blank, then we want to also make the targets...
# CHILDREN:   NON-BLANK
# targets :   NON-BLANK
# Note that it makes the children FIRST, then the local files. So it makes the most deep files first, then the shallower ones
all: all_children
	$(MAKE) $(targets); \

try:
	$(MAKE) $(early_targets)
	$(MAKE) all_children
	$(MAKE) $(targets)

clean: clean_children
	rm -f $(files_to_clean) $(wildcard *.tmp); \

back: back_children
	mkdir -p Backup; \
	cp -f $(things_to_back_up) $(wildcard Makefile*) Backup; \

else   # NESTED IF 2: else (if targets IS IN FACT blank)
# If we get here, then CHILDREN has some entries, but targets is blank...
# CHILDREN:   NON-BLANK
# targets :   BLANK
all: all_children
clean: clean_children
back: back_children
early_targets: early_targets
children: all_children
subdirs: all_subdirs

endif  # END OF NESTED IF 2 (done checking to see if targets is blank or not)
# If we get here, CHILDREN has something in it, and targets can be anything
# CHILDREN:   NON-BLANK
# targets :   BLANK  or  NON-BLANK

# Note: the old make rule won't handle an enormous number of commands ("argument list too long" error)
# Hence this new rule.
make:
	for dir in $(CHILDREN); do \
		mkdir -p "$$dir" ;\
		cd "$$dir" ;\
		rm -f Makefile ;\
		ln -s $(CHILD_MAKEFILE) Makefile ;\
		$(MAKE) make ;\
		cd .. ;\
	done


remote: remote_children
	mkdir -p $(REMOTE_DIR); \
	$(foreach u, $(URL), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) '$(u)'; \
	) \

remote_clean: remote_clean_children
	rm -rf $(REMOTE_DIR); \

back_clean: back_clean_children
	rm -rf Backup; \


# This lets you apply ANY rule recursively to the subdirectories.
# It is extremely convenient!
# If you have a rule named "plot" in every sub-directory (each of the CHILDREN),
# you can call it by saying "make plot_children"
%_children:
	for dir in $(CHILDREN); do \
		mkdir -p "$$dir"  ;\
		cd "$$dir"    ;\
		$(MAKE) $*  ;\
		cd .. ;\
	done

# This appears to be an alias for some unknown reason.
%_subdirs: %_children

else   # NESTED IF 1
# (CHILDREN is empty if we get here)
# CHILDREN:   BLANK
# targets :   BLANK  or  NON-BLANK

# If we get here, CHILDREN is empty, but we might still have targets
all: $(targets)

clean:
	rm -f $(files_to_clean) $(wildcard *.tmp); \

back:
	mkdir -p Backup; \
	cp -f $(things_to_back_up) $(wildcard Makefile*) Backup; \

back_clean:
	rm -rf Backup; \

make:

remote:
	mkdir -p $(REMOTE_DIR);
	$(foreach u, $(URL), \
	   $(REMOTE_FTP) -P $(REMOTE_DIR) '$(u)'; \
	)

# end of empty CHILDREN
endif # END OF NESTED IF 1 (done checking the status of CHILDREN)


# Anything below here is always the same no matter what CHILDREN and targets are
what_settings:
	echo 'CHILDREN is ' $(CHILDREN)
	echo 'early_targets is ' $(early_targets)
	echo 'targets is ' $(targets)
	echo 'side_effects is ' $(side_effects)
	echo 'CHILD_MAKEFILE is ' $(CHILD_MAKEFILE)

