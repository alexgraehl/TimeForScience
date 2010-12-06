# This makefile for filtering link data by cluster support
# It should live in the dataset's ClusterFilter directory.

MODULES_FILE = ../lists_t.tab

LINKS_FILE = ../../data.tab.gz

TARGETS = links.tab.gz

clean:
	rm ${TARGETS}

links.tab.gz:
	cluster_supported_links.pl -d -l ${LINKS_FILE} -c ${MODULES_FILE} \
	| cap.pl \#Element1,Element2,Correlation,Dimensions,Fisher-Z,Rank,ValidNeighbors,clusters\
	> 1.tmp ;
	\
	mv 1.tmp links.tab ;
	\
	gzip links.tab ;
	