default: all

# READ THIS IF YOU ARE HAVING PROBLEMS:
# ** Remember that the definition of CHILDREN/targets/etc
# ** should typically come BEFORE recurse.mak is included!
#
# i.e.
# targets = file.tab
# include recurse.mak   <-- (this will work)

ifeq ($(origin targets), undefined)
ifeq ($(origin CHILDREN), undefined)
ifeq ($(origin CHILD_MAKEFILE), undefined)
ifeq ($(origin early_targets), undefined)
$(error "********** MAKEFILE PROBLEM ************* You need to define at least one of these variables: targets, CHILDREN, or early_targets  when you include recurse.mak. Probably the problem is that you did not define these variables BEFORE you included recurse.mak! You must define the variables first, and THEN include recurse.mak!!! (It actually usually works to define them later, but any ifeq statements in recurse.mak that rely on them will break.)\n")
endif
endif
endif
endif


# This is a clone of quick.mak by Alex Williams. It has some additional features and error checking.
# **************
# Overview:
# This is a Makefile that lets you recursively call other makefiles in subdirectories without having to figure it out yourself every single time.
# **************

# **************
# USAGE:
# You can often just replace "include quick.mak" with "include recurse.mak" , and since the same variable names are used, it should just work.
# Basically, you just set some variables (as described below) and then include recurse.mak.
# **************

# **************
# Additional Error Checking:
# When recursing, it checks to see if files exist or not (in the RUN_RECURSIVELY_ONLY_IN_THESE_DIRS perl code).
# (quick.mak has no way of checking to see if Makefiles exist or not).
# **************

# **************
# The order of operations for recurse.mak:
#  early_targets: These files (located in the current directory) get made FIRST. Maybe the SUBDIRS/CHILDREN depend on these files.
#  CHILDREN: Then these subdirectories are traversed
#  SUBDIRS: Then *these* subdirectories are traversed
#  targets: Then finally these files (located in the current directory) are made. These might depend on the SUBDIRS/CHILDREN.
# **************

# **************
# Variables that you'll need to set in the makefile that includes recurse.mak:
# NOTE: Set these BEFORE (i.e., earlier in the file) you include recurse.mak!!!
# early_targets   <-- files to be made *before* recursing to CHILDREN/SUBDIRS
# targets         <-- files to be made *after* traversing CHILDREN/SUBDIRS
# CHILDREN        <-- these directories are created if you run "make make" 
#                     (manually-generated / important files should NEVER EVER be
#                     in a CHILDREN directory--they could be auto-deleted by make_clean)
# CHILD_MAKEFILE  <-- this is the makefile that is linked in each CHILDREN directory if you run "make make"
# SUBDIRS (rare)  <-- these are directories that are traversed, but shouldn't be deleted.
#                     A CHILD_MAKEFILE will be linked here if there is not already a Makefile in the SUBDIR
# **************

# **************
# How to run it:
#
# If you have defined any CHILDREN, you should first type:   make make
#  in order to set up the directory structure defined by CHILDREN and CHILD_MAKEFILE
#
# To run the rules, type:    make    in order to make all the targets (and targets for child makefiles).
# **************

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
# Note that we process the CHILDREN directories *BEFORE* we handle the targets.
# Important files and real makefiles should NEVER EVER be in a CHILDREN directory, as they are subject to automatic
# deletion. Only generated files and symlinks to important files should be in here.
CHILDREN ?=

# SUBDIRS is a list of subdirectories that we call and run commands in, but they are NOT ever cleaned.
# Although we will link the CHILD_MAKEFILE in these if no makefile exists, SUBDIRS are really intended
# to be under manual control, rather than the automatic generation / deletion of CHILDREN.
SUBDIRS ?=
BACKUP_FOLDER ?= Backup
REMOTE_FOLDER ?= Remote

# "early_targets" are the files we generate BEFORE running the children
early_targets ?=

# Targets are files that we could generate with a "make (this target)" command. Note that these are generated AFTER we handle the CHILDREN
targets ?=
# side_effects are files that get generated when we're making the targets. We keep track of them here so we'll know to delete them when a user runs "make clean". But we never "make (side_effect)" any of them.
side_effects ?=

files_to_clean    = $(early_targets) $(targets) $(side_effects)
things_to_back_up = $(early_targets) $(targets)

# Below: Set the default CHILD_MAKEFILE in case none is defined.
# Note that the path is from the perspective of the subdirectory (i.e., everything requires one more level of ../ than you might expect)
CHILD_MAKEFILE ?= ../Lib/child.mak

# "PHONY" targets are all rules that aren't also the name of a file.
# e.g. the "all" rule does not generate a file named "all," so it's a phony target
# (if you don't specify these, then a file named "all" can sabotage "make all" from running.
# You will get an error like "make: "all" is already up-to-date")
.PHONY: default all local eradicate make_clean clean back remote remote_clean back_clean make vars

URL     ?=
url_ndx  = $(shell seq 1 $(words $(URL)))
# The URL-related information is for automatically downloading "Remote" data to our local drive, then using it in some makefile step.


# Make the early_targets first, THEN the children, THEN finally the targets
all: vars
	@perl -e 'if (qq{} !~ qq{$(early_targets)}) { system(qq{$(MAKE) $(early_targets)}); }'
	$(call RUN_RECURSIVELY,all)
	@perl -e 'if (qq{} !~ qq{$(targets)}) { system(qq{$(MAKE) $(targets)}); }'

# Just make the targets in THIS directory, don't recursively do the children. Sometimes convenient for updating summary files without having to check all the prerequisites.
local:
	$(MAKE) $(early_targets) $(targets)

# Argument 1: what rule do we want to make in each directory?
RUN_RECURSIVELY = \
	$(call RUN_RECURSIVELY_ONLY_IN_THESE_DIRS,$(1),$(CHILDREN) $(SUBDIRS),0)

# Set RECURSE_STOP_ON_WARNINGS to 1 in order to stop when we are missing a makefile or child directory.
# Set it to 0 to continue on anyway.
RECURSE_STOP_ON_WARNINGS ?= 1

# Argument $(1): what rule do we want to make in each directory?
# Argument $(2): which directories do we want to look in?
# Argument $(3): whether ot not to create nonexistent directories and link the child makefiles (used for "make make"). 1 = create them, 0 = don't.
#  If $(3) is set to any true value, then we create directories and link the child makefile for each one, and then run whatever rule was being called recursively.
#  This only really makes sense for "make make," which is the rule that sets up the CHILDREN/SUBDIRS directory structure.
# NOTE: The @perl means "do not echo this rule to the command line when this rule is run"
RUN_RECURSIVELY_ONLY_IN_THESE_DIRS = \
	@perl -e '\
	use strict; use warnings; \
	my $$shouldMakeChildrenDirs = $(3); \
	my $$allDirStr = "$(2)"; \
	my @allDirArray = split(/\s+/,$$allDirStr); \
	foreach my $$dir (@allDirArray) { \
		if ($$shouldMakeChildrenDirs) { \
			print STDERR qq{recurse.mak: Making directory $$dir...\n}; \
			system(qq{mkdir -p $$dir}); \
		} \
		my $$prevDir = `pwd`; chomp($$prevDir); \
		if (length($$dir) == 0) { next; } \
		if (not -d qq{$$dir}) { \
			print STDERR qq{recurse.mak: WARNING: No directory <$$dir> was found (although this entry was in the CHILDREN / SUBDIRS for this Makefile).\n}; \
			if ($(RECURSE_STOP_ON_WARNINGS)) { die "recurse.mak: Quitting due to this error. To not stop on errors, set \"RECURSE_STOP_ON_WARNINGS = 0\" before including recurse.mak.\n"; } \
			next; \
		} \
		print STDERR qq{recurse.mak: Recursively running <$(MAKE) $(1)> from directory <$$prevDir>.\n}; \
		print STDERR qq{recurse.mak: Changing into directory <$$dir>\n}; \
		chdir qq{$$dir}; { \
			my $$thisDir = `pwd`; chomp($$thisDir); \
			\
			if ($$shouldMakeChildrenDirs) { \
				if ((not -f qq{Makefile}) or (-z qq{Makefile})) { \
					print STDERR qq{recurse.mak: Linking the child makefile <$(CHILD_MAKEFILE)> to the file <$$thisDir/Makefile>.\n}; \
					system(qq{ln -f -s $(CHILD_MAKEFILE) Makefile}); \
					if ((not -f qq{$(CHILD_MAKEFILE)})) { \
						print STDERR qq{recurse.mak: WARNING: Sym-linked <$$thisDir/Makefile> to the file <$(CHILD_MAKEFILE)>--but this file does not (yet) exist.\n}; \
						if ($(RECURSE_STOP_ON_WARNINGS)) { die "\n\n*************************\nrecurse.mak: Quitting due to the missing Makefile. You should make sure a make-compatible file exists at <$(CHILD_MAKEFILE)>. To not stop on errors, set the variable \"RECURSE_STOP_ON_WARNINGS = 0\" in the Makefile that includes recurse.mak.\n*************************\n\n"; } \
					} \
				}	\
			} \
			\
			if (-f qq{Makefile} && not (-z qq{Makefile})) { \
				print STDERR qq{recurse.mak: Recursing with rule <$(1)> in <$$thisDir/Makefile>...\n}; \
				my $$returnCode = system(qq{$(MAKE) $(1)}); \
				if ($$returnCode != 0) { \
					die qq{recurse.mak: ERROR: <$(MAKE) $(1)> exited with non-zero exit code <$$returnCode>. Quitting recursive make.\n}; \
				} \
			} else { \
				print STDERR qq{recurse.mak: WARNING: Skipping <$$thisDir>, which has no Makefile.\n}; \
				if ($(RECURSE_STOP_ON_WARNINGS)) { die "\n\n*************************\nrecurse.mak: Quitting due to the missing Makefile. To not stop on errors, set the variable \"RECURSE_STOP_ON_WARNINGS = 0\" in the Makefile that includes recurse.mak.\n*************************\n\n"; } \
			} \
		} chdir($$prevDir); \
	} \
	'


PERLCLEAN = perl -e "if (-d q{$$dir}) { \
					system(qq{/bin/rmdir  --ignore-fail-on-non-empty } . q{\"$$dir\"} ); \
				} \
				else { \
					system(qq{/bin/rm -f } . q{\"$$dir\"} ); \
				} \
				"

clean:
	@for dir in $(files_to_clean) $(wildcard *.tmp) $(wildcard *.temp); do \
		$(PERLCLEAN) ; \
	done ;

cleaner:
	$(MAKE) clean ;
	$(call RUN_RECURSIVELY,cleaner) ;

# COMPLETELY REMOVE all the child directories. Dangerous! Note: should NOT remove $(SUBDIRS). Children are auto-generated, SUBDIRS are not.
make_clean cleanest: # This was a rule in quick.mak. COMPLETELY REMOVES all the child directories. Note that it does NOT remove $(SUBDIRS). Can be dangerous!
	$(RM) -r $(CHILDREN) ;
	$(RM) -r $(files_to_clean) ;
	$(RM) $(wildcard *.tmp) $(wildcard *.temp) ;


back:
	$(call RUN_RECURSIVELY,back)
	mkdir -p $(BACKUP_FOLDER) ;
	cp $(things_to_back_up) $(wildcard Makefile*) $(wildcard *.mak)    $(BACKUP_FOLDER) ;

remote:
	$(call RUN_RECURSIVELY,remote)
	mkdir -p $(REMOTE_FOLDER); \
	$(foreach u, $(URL), \
	   $(REMOTE_FTP) -P $(REMOTE_FOLDER) '$(u)'; \
	) \

remote_clean:
	$(RM) -rf $(REMOTE_FOLDER); \

back_clean:
	$(RM) -rf $(BACKUP_FOLDER); \

# Generates the CHILDREN directories and then runs MAKE. Same as "make make ; make"
gen generate:
	$(MAKE) make ;
	$(MAKE) ;

# For the rule "make": Note that "cd -" is in place because "cd .." only moves up ONE level, whereas it is possible that the CHILDREN is actually a chain of *more than one* directory.
# In quick.mak, it was originally "cd .." .
make:
	$(call RUN_RECURSIVELY_ONLY_IN_THESE_DIRS,$@,$(CHILDREN) $(SUBDIRS),1)

old_recurse:
	for dir in $(CHILDREN) $(SUBDIRS); do \
		mkdir -p "$$dir"		;\
		cd "$$dir"			;\
		  ln -s $(CHILD_MAKEFILE) Makefile	;\
		  $(MAKE) $@				;\
		cd -				;\
	done ;


# Anything below here is always the same no matter what CHILDREN and targets are
vars:
	@echo "-------------------------------------------------"
	@echo "Recurse.mak is running..."
	@echo '   early_targets: ' $(early_targets)
	@echo '         targets: ' $(targets)
	@echo '        CHILDREN: ' $(CHILDREN)
	@echo '    side_effects: ' $(side_effects)
	@echo '  CHILD_MAKEFILE: ' $(CHILD_MAKEFILE)

#.cvsignore:
#	echo $(files_to_clean) | transpose.pl -d ' ' > $@;

