
# Below: "guard" statement: Don't include this file multiple times, or we get annoying warnings
ifndef FREQUENT_FUNCTIONS_MAKEFILE_INCLUDED
FREQUENT_FUNCTIONS_MAKEFILE_INCLUDED = 1

# This file consists of a list of really common things that we often want to do in makefiles.

# Print the innermost directory
# TAB is a literal tab
EMPTY =
TAB   = $(EMPTY)	$(EMPTY)
THIS_DIR = $(shell pwd | sed 's/.*\///')

#perl -e 'my @x = split("/",`pwd`); print $x[-1]'
INNERMOST_DIR = $(THIS_DIR)

F_REMOVE_HEADER_LINE = tail -n +2
f_REMOVE_HEADER_LINE = ${F_REMOVE_HEADER_LINE}
FUNC_REMOVE_HEADER_LINE = $(F_REMOVE_HEADER_LINE)



F_REMOVE_TRAILING_WHITESPACE = sed 's/[ 	]\+$$//'
FUNC_REMOVE_TRAILING_WHITESPACE = $(F_REMOVE_TRAILING_WHITESPACE)

F_REMOVE_LINES_WITH_ONLY_WHITESPACE = grep -v '^[ 	]*$$'
FUNC_REMOVE_LINES_WITH_ONLY_WHITESPACE = $(F_REMOVE_LINES_WITH_ONLY_WHITESPACE)

F_REMOVE_LEADING_AND_TRAILING_WHITESPACE = sed 's/^[ 	]*//;s/[ 	]*$$//'

# Run this like so:   cat MYFILE | $(call F_PREPEND_LINE_NUMBERS) > OUTPUTFILE
# Equivalent to "lin.pl" for cases without special arguments to lin.pl
F_PREPEND_LINE_NUMBERS = awk '{ print FNR "\t" $$0 }'


# Below: This works a lot like fasta2single_line.pl (but is comma-safe... safer, anyway)
# Note: sometimes you get modules with NO results, in which case you just have two
# consecutive lines like
# >Module
# >Module2
# That's why at the end sometimes you have blank entries
FUNC_CONVERT_SETS_OVERLAP_OUTPUT_TO_TAB_DELIMITED = \
		perl -p -e 's/^>/==LINEBREAK==/g' 	\
		| perl -p -e 's/\n/\t/g'			\
		| perl -p -e 's/==LINEBREAK==/\n/g'	\
		| flatten.pl					\
		| sed 's/\\,/==COMMA==/g'		\
		| perl -p -e 's/\t([^,]+),([^,]+),/\t$$1\t$$2\t/g' \
		| sed 's/==COMMA==/\\,/g' \
		| select.pl -k 2 -lne ''


FUNC_SWITCH_SET_MAJOR_WITH_ELEMENT_MAJOR = \
	flatten.pl \
		| expand.pl -k 2 \
		| sort -k 1,1 -t '	'

# Replaces:
#   tab  (nothing here)  tab
# With:
#   NA   tab       NA         tab    NA
FUNC_REPLACE_EMPTY_WITH_NA = \
	sed -e 's/\t\t/\tNA\t/g' -e 's/\t\t/\tNA\t/g' -e 's/\t$$/\tNA/g' -e 's/^\t/NA\t/g'


# Prints lines EXCEPT lines where column 1 is equal to column 2
# $(1): number/index of first column (counting from 1 = first column... not [0]!)
# $(2): number/index of second column (counting from 1 = first column... not [0]!)
# $(3): The thing to print for each skipped row (you can remove these later with grep)
# Example usage: $(call NO_SELF_COMPARISONS,1,2,SELF_COMPARISON)
NO_SELF_COMPARISONS = \
	perl -n -e 'chomp($$_); my @x = split("\t",$$_); if (scalar(@x) < ($(2))) { print "$$_\n"; } else { if ($$x[$(1)-1] ne $$x[$(2)-1]) { print "$$_\n"; } else { print "$(3)\n"; } }'

# Prepends a count of how many tab-delimited items are on a line, for each line.
# $(1): input   $(2): output
PREPEND_COUNT_TO_LIST_FILE = \
	cat $(1) \
		| row_stats.pl -h 0 -count \
		| cut -f 2 \
		| tail -n +2 \
		| paste - $(1) \
			> $(2) ;




# 1: thing to assign (does not actually assign it--you call this function, and it uses the old value if you didn't assign it to anything new)
# 2: new value if we are in this subdir
# 3: the subdir
#FUNC_ASSIGN_IF_WE_ARE_IN_THIS_SUBDIR = \
#	$(shell perl -e 'use strict; use warnings; \
#		if (`pwd` =~ /\/$(3)\//) { \
#			print q{$(2)}; \
#		} else print q{$(1)};




endif  # end of "guard" at top of file
