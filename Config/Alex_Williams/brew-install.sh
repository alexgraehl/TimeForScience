#!/bin/bash


<<<<<<< HEAD
GAMES="$CKAN_KERBAL"
BIOINFORMATICS="fastqc tophat bowtie2 bedtools samtools htslib bcftools vcflib vcftools"
=======
brew tap homebrew/science
#brew tap homebrew/dupes
# =================================
>>>>>>> b0d4bb53c3491419186dca4e271cbd4df270e109

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


