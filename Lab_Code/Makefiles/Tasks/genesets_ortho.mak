CHILDREN       = $(shell grep '^ORGANISMS\>' ~/Map/Makefile.common | cut -f 2 -d =)
CHILD_MAKEFILE = $(HOME)/Map/Templates/Make/genesets_ortho_child.mak

include $(MAPDIR)/Templates/Make/quick.mak


