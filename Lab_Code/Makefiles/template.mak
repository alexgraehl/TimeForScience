SHELL     := /bin/bash
EMPTY     :=
SPACE     :=$(EMPTY) $(EMPTY)
TAB       :=$(EMPTY)	$(EMPTY)
# Above: TAB: A literal tab! Make sure not to change the whitespace on this line (the tab is, naturally, invisible)
# IMMEDIATE EVALUATION>>  := means "set this variable to be whatever the value on the right side is at this very moment"
# LAZY EVALUATION>>   =  means "set this variable to the equation on the right side, and then evaluate the right side whenever it comes up.

downloaded=
generated=test_1


all: $(generated) $(downloaded)

clean:
	${RM} $(generated)

eradicate:
	${RM} $(generated)
	${RM} $(downloaded)

# ====================== REAL RULES BELOW HERE =============

test_1:
	@echo "I should be a rule"
