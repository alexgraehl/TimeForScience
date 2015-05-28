#!/bin/bash -
#
# To use this filter with less, define LESSOPEN:
# export LESSOPEN="|/usr/bin/lesspipe.sh %s"
# (Note that that is the default one, not this advanced one)

# In order to use this fancy less script with the command "vv" (while keeping less
# unchanged), add the following line to your .cshrc file. Bash users will have to
# figure out how to make an alias on their own systems.
# alias vv "env LESSOPEN='|${MYSRC}/lab_apps/shell_scripts/lesspipe_advanced.sh %s' less -S --LINE-NUMBERS --status-column --ignore-case -R -f \!*"

# Changes from the default script are annotated below

smartdecompress() {
    # "Smart" picking of a decompressing utility (or none!) and echoing to STDOUT. Normally you pipe this into ANOTHER utility!
    # Example:   smartdecompress "$1" | sed 's/ABC/DEF/'
    case "$1" in
	*.gz)	DECOMPRESSOR="gunzip -c" ;;
	*.bz2)	DECOMPRESSOR="bunzip2 -c" ;;
	*)	DECOMPRESSOR="cat" ;;
    esac
    $DECOMPRESSOR $1; # <-- this is what we pass into the relevant "smart" function below
}

#STO="\033[1m" # start color... not useful?
#STA="\033[31m" # start color A: 31 = red
#STC="\033[36m" # start color C: 36 = cyan. 32 = green. 35 = magenta
#STG="\033[33m" # start color G: 33 = yellow
#STT="\033[35m" # start color T: 36 = cyan. 35 = magenta
#STN="\033[45m" # start color N: magenta background, doesn't change the foreground
#RES="\033[0m" # RESET color
basecolor() {
    ## Colors the A, C, G, and T
    ## Assumes that you pass something INTO it via a pipe -- otherwise it fails
    ## Example:  cat myfile | sed 's/1/2/' | basecolor 
    #sed -e 's/\(AA*\)/\1'"$(printf ${STA}A${RES})"'/g' -e 's/C/'"$(printf ${STC}C${RES})"'/g' -e 's/G/'"$(printf ${STG}G${RES})"'/g' -e 's/T/'"$(printf ${STT}T${RES})"'/g' -e 's/N/'"$(printf ${STN}N${RES})"'/g'
    perl ${TIME_FOR_SCIENCE_DIR}/Config/1_Shell_Config/lesspipe_extra_bamfile_coloration.pl
    #sed -e 's/A/'"$(printf ${STA}A${RES})"'/g' -e 's/C/'"$(printf ${STC}C${RES})"'/g' -e 's/G/'"$(printf ${STG}G${RES})"'/g' -e 's/T/'"$(printf ${STT}T${RES})"'/g' -e 's/N/'"$(printf ${STN}N${RES})"'/g'

    # Colorize based on entire blocks of identical characters. Turns out not to be any faster, unfortunately!
    #sed -e 's/\(AA*\)/'"$(printf $STA)"'\1'"$(printf $RES)"'/g' -e 's/\(CC*\)/'"$(printf $STC)"'\1'"$(printf $RES)"'/g' -e 's/\(GG*\)/'"$(printf $STG)"'\1'"$(printf $RES)"'/g' -e 's/\(TT*\)/'"$(printf $STT)"'\1'"$(printf $RES)"'/g' -e 's/\(NN*\)/'"$(printf $STN)"'\1'"$(printf $RES)"'/g'
}

lesspipe() {
    case "$1" in
	*.[1-9n]|*.man|*.[1-9n].bz2|*.man.bz2|*.[1-9].gz|*.[1-9]x.gz|*.[1-9].man.gz) ## <--- For viewing multi-part archives!
	    case "$1" in
		*.gz)	DECOMPRESSOR="gzip  -dc" ;;
		*.bz2)	DECOMPRESSOR="bzip2 -dc" ;;
		*)	DECOMPRESSOR="cat" ;;
	    esac
	    if $DECOMPRESSOR -- "$1" | file - | grep -q troff; then
		if echo "$1" | grep -q ^/; then	# absolute path
		    man -- "$1" | cat -s
		else
		    man -- "./$1" | cat -s
		fi
	    else
		$DECOMPRESSOR -- "$1"
	    fi ;; ## <-- Done with the "multi-part archive" part
  # Alex: CHANGED from default lesspipe: Now we load any files ending in
  # .tab or .matrix (including gzipped files) and run them through
  # "sheet.pl", which will format them as a quasi-spreadsheet with
  # columns aligned.   .csv files could be added in a similar fashion
  # by just changing the delimiter used by sheet.pl
  # For a faster but less feature-filled
  #  UNIX-tool-only version of this, instead of "sheet.pl", use "column -t -s '	'"
	*.tab.gz|*.matrix.gz) gunzip -c "$1" | sheet.pl --notify --color="always" ;;
	*.tab|*.matrix) sheet.pl --notify --color="always" "$1" ;; # Maybe should use sheet.py instead?
	
	*.bam) if [[ -x "`which samtools`" ]]; then  # Color code the bases A, C, G, and T (and N)
		   echo -e '##################\n# Viewing a BAM file with "samtools view -h" -- The actual file is in a binary format #\n##################' ; samtools view -h "$1" | basecolor ## Use Samtools to view a BAM ("binary sam") RNA-seq file
	       else  ## Use Samtools to view a BAM ("binary sam") RNA-sequence file
		   echo -e "ERROR: LESSPIPE_ADVANCED CANNOT VIEW THIS FILE, BECAUSE 'samtools' IS NOT INSTALLED or CANNOT BE FOUND.\nPlease install samtools (or make sure it is on your path and can be executed) and try again!"
	       fi ;;
	###########	*.bam) echo '######\n### Viewing a BAM file with "samtools view -h" -- The actual file is in a binary format ###\n######' ; samtools view -h "$1" | sed 's/A/'"$(printf '\033[1mA\033[0m')"'/g' ;; ## Use Samtools to view a BAM ("binary sam") RNA-sequence file
	*.sam|*.sam.gz|*.sam.bz2|*.dna.txt|*.rna.txt|*.dna.txt.gz|*.rna.txt.gz|*.seq.gz|*.seq.txt|*.seq.txt.gz|*.seq|*.fasta|*.fasta.gz|*.fasta.bz2|*.fastq|*.fastq.gz|*.fastq.bz2|*.fa|*.fa.gz|*.fa.bz2|*.fq|*.fq.gz|*.fq.bz2) smartdecompress "$1" | basecolor ;; ## Use Samtools to view a BAM ("binary sam") RNA-sequence file
	*.tar) tar tvvf "$1" ;;
	*.tgz|*.tar.gz|*.tar.[zZ]) tar tzvvf "$1" ;;
	*.tar.bz2|*.tbz2) bzip2 -dc "$1" | tar tvvf - ;;
	*.[zZ]|*.gz) gzip -dc -- "$1" ;;
	*.bz2) bzip2 -dc -- "$1" ;;
	*.zip) zipinfo -- "$1" ;;
	*.rpm) rpm -qpivl --changelog -- "$1" ;;
	*.cpi|*.cpio) cpio -itv < "$1" ;;
	*.gif|*.jpeg|*.jpg|*.pcd|*.png|*.tga|*.tiff|*.tif)
	    if [ -x "`which identify`" ]; then
		identify "$1"
	    else
		echo "No identify available"
		echo "Install ImageMagick to browse images"
	    fi ;;
	*)
    esac
}

if [ -d "$1" ] ; then
    /bin/ls -alF -- "$1"
else
    lesspipe "$1" 2> /dev/null
fi
