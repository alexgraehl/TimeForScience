#!/bin/sh -
#
# To use this filter with less, define LESSOPEN:
# export LESSOPEN="|/usr/bin/lesspipe.sh %s"
# (Note that that is the default one, not this advanced one)

# In order to use this fancy less script with the command "vv" (while keeping less
# unchanged), add the following line to your .cshrc file. Bash users will have to
# figure out how to make an alias on their own systems.
# alias vv "env LESSOPEN='|${MYSRC}/lab_apps/shell_scripts/lesspipe_advanced.sh %s' less -S --LINE-NUMBERS --status-column --ignore-case -R -f \!*"

# Changes from the default script are annotated below

lesspipe() {
  case "$1" in
  *.[1-9n]|*.man|*.[1-9n].bz2|*.man.bz2|*.[1-9].gz|*.[1-9]x.gz|*.[1-9].man.gz)
	case "$1" in
		*.gz)	DECOMPRESSOR="gunzip -c" ;;
		*.bz2)	DECOMPRESSOR="bunzip2 -c" ;;
		*)	DECOMPRESSOR="cat" ;;
	esac
	if $DECOMPRESSOR -- "$1" | file - | grep -q troff; then
		if echo "$1" | grep -q ^/; then	#absolute path
			man -- "$1" | cat -s
		else
			man -- "./$1" | cat -s
		fi
	else
		$DECOMPRESSOR -- "$1"
	fi ;;

  # Alex:
  # CHANGED from default lesspipe: Now we load any files ending in
  # .tab or .matrix (including gzipped files) and run them through
  # "sheet.pl", which will format them as a quasi-spreadsheet with
  # columns aligned.   .csv files could be added in a similar fashion
  # by just changing the delimiter used by sheet.pl
  # For a faster but less feature-filled
  #  UNIX-tool-only version of this, instead of "sheet.pl", use "column -t -s '	'"
  *.tab.gz|*.matrix.gz) gunzip -c "$1" | sheet.pl --notify ;;
  *.tab|*.matrix) sheet.pl --notify "$1" ;;

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
