CHILDREN       = $(shell grep '^ORGANISMS\>' $(MAPDIR)/Makefile.common | cut -f 2 -d =)
CHILD_MAKEFILE = $(HOME)/Map/Templates/Make/interactions_ortho_child.mak

include $(MAPDIR)/Templates/Make/quick.mak


