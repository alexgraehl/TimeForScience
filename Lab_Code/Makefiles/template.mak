SHELL     := /bin/bash
EMPTY     :=
SPACE     :=$(EMPTY) $(EMPTY)
TAB       :=$(EMPTY)	$(EMPTY)

# Above: TAB: A literal tab! Make sure not to change the whitespace on this line (the tab is, naturally, invisible)
# IMMEDIATE EVALUATION>>  := means "set this variable to be whatever the value on the right side is at this very moment"
# LAZY EVALUATION>>   =  means "set this variable to the equation on the right side, and then evaluate the right side whenever it comes up.

# .ONESHELL:
# .SHELLFLAGS := -eu -o pipefail -c
# .DELETE_ON_ERROR:
# MAKEFLAGS += --warn-undefined-variables
# MAKEFLAGS += --no-builtin-rules


irreplacable=
downloaded=
generated=test_1

all: $(generated) $(downloaded)

clean:
	@echo "Deleting the Makefile-generated files."
	${RM} $(generated)

eradicate:
	@echo "Deleting the Makefile-generated files."
	${RM} $(generated)
	@echo "Additionally, deleting the downloaded files. Hopefully they can be downloaded again, but no guarantees!"
	${RM} $(downloaded)
	@echo "Not deleting anything from the 'irreplacable' list"

# ====================== REAL RULES BELOW HERE =============

test_1:
	@echo "I should be a rule"
