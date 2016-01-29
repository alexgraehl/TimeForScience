#!/bin/bash


brew tap homebrew/science
#brew tap homebrew/dupes
# =================================

# =================================
# Less common bioinformatics tools
CLUSTAL="homebrew/science/clustal-w" # clustalw2
MIN_ION="poretools"
AMAZON="awscli"
BLAST="blast"
PHYLIP_TOOLS="phylip"
# =================================
GPG_TOOLS="Caskroom/cask/gpgtools"
# =================================
UTILITIES="homebrew/dupes/unzip"
# =================================
PROGRAMMING="bash wget emacs ess git mercurial"
# =================================
BIOINFORMATICS="fastqc tophat bowtie2 bedtools samtools bcftools"
# =================================
brew tap homebrew/games
CKAN_KERBAL="ckan"   # Kerbal space program package manager
GAMES="$CKAN_KERBAL"
# =================================

brew install $UTILITIES $PROGRAMMING $BIOINFORMATICS $GAMES


