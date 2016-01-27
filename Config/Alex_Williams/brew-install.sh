#!/bin/bash

brew tap homebrew/games
brew tap homebrew/science


CKAN_KERBAL="ckan"

GAMES="$CKAN_KERBAL"
BIOINFORMATICS="fastqc tophat bowtie2 bedtools samtools"

MIN_ION="poretools"
AMAZON="awscli"
BLAST="blast"
PHYLIP_TOOLS="phylip"
GPG_TOOLS="Caskroom/cask/gpgtools"
CLUSTAL="homebrew/science/clustal-w" # clustalw2

PROGRAMMING="bash wget emacs ess git mercurial"

brew install $GAMES $PROGRAMMING $BIOINFORMATICS
