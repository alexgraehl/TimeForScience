#
#

targets = et2.txt
intermediates = 
keep_us_please = etym_parsed_1.txt enwiktionary-latest-pages-articles.xml

all: et2.txt
	echo "Done?"


enwiktionary-latest-pages-articles.xml: 
	echo "If this file is not around, you will need to download it. It is 200 MB when compressed."
	echo "You need to download the file named: enwiktionary-latest-pages-articles.xml"
	echo " from this web site: http://dumps.wikimedia.org/enwiktionary/latest/"
	echo "It is 200 MB when uncompressed. Then simply bunzip2 <filename> to get this xml file."


#~/TimeForScience/Self-contained_Projects/Alex_Williams/Etymology_Plugin/alex_parse_enwiktionary-latest-pages-articles.pl:
#

etym_parsed_1.txt: enwiktionary-latest-pages-articles.xml
	perl ~/TimeForScience/Self-contained_Projects/Alex_Williams/Etymology_Plugin/alex_parse_enwiktionary-latest-pages-articles.pl $< > $@

## Just take out columns 2 and 3, those are the only ones we want.
## Column 1 is just the line number!
et2.txt: etym_parsed_1.txt
	cut -f 2,3 $< > $@
	
## Strip out everything that isn't a wWORD or an etymology (etyl)
et3.txt: et2.txt
	cat $< | perl -e '\
while(my $$line = <>) {	\
	chomp($$line);		\
	if ($$line =~ m/WORD\S*\s+(\S+)/) {	\
		print "\n" . $$1 . "\t";	\
	}	\
	my @a = ($$line =~ m/({{etyl.*?}})/g);	\
	if(scalar(@a) > 0) {	\
		print join("\t", @a);	\
		print "\n";	\
	} \
} \
' > $@


## perl multi-line replace: search for '{{' starting a line,
## and then put that back on the previous line (by removing the newline!)
et4.txt: et3.txt
	perl -00pe 's/\n[{][{]/\t{{/gxms' $< > $@

## Take out anything with no listed etymology.
## grep -P means PERL COMPATIBLE MODE (much faster than regular grep)
## Look for lines that have {{etyl on them ONLY. Discard the others.
et5.txt: et4.txt
	grep --perl-regexp '{{etyl' $< > $@


## Find entries where the word has AT LEAST ONE number or
## Latin letter in it. Removing all words in 100% foreign
## character sets.
## Regular expression explanation:
## starts with anything, then some letter we recognize, then any set of non-whitespace letters, then a tab
et6.txt: et5.txt
	grep --perl-regexp '^\S*[a-zA-Z0-9]\S*\t' $< > $@

## export LC_ALL='C' works around the fact that sort
## has issues dealing with certain foreign byte sequences when
## you are doing --ignore-case sorts.
et7.txt: LC_ALL = 'C'
et7.txt: et6.txt
	OLDLC=$$LC_ALL
	export LC_ALL='C' && sort --ignore-case $< > $@ ## no idea why, but EXPORT has to be on this same line as the sort!
	export LC_ALL=$$OLDLC

## Ok, take out anything that starts with a hyphen, quote, or &.
## or ENDS with a hyphen
et8.txt: et7.txt
	cat $< | grep --perl-regexp -v "^[-&']" \
		| grep --perl-regexp -v "[-]$$" \
			> $@

## Use Josh's script to collapse duplicate entries!
et9.txt: et8.txt
	cat $< | expand.pl > $@

## Finally, remove the "etyl" and brackets
eta.txt: et9.txt
	cat $< | sed -e 's/{{etyl|//g' -e 's/}}//g' > $@


j1.txt: eta.txt
	echo 'const WORD = {' > $@
	perl -n -e '\
chomp($$_);	\
my @x = split(/\t/, $$_);	\
print qq{"$$x[0]" : [} . qq{"$$x[1]"}. qq{, "COLORPLACEHOLDER"],\n};	\
' $< >> $@
	echo '};' >> $@

j2.txt: j1.txt
	cat $< \
		| perl -p -e 's/(.* : [[]\".*en.*?\".*)COLORPLACEHOLDER/$${1}#DDDDFF/' \
		| perl -p -e 's/(.* : [[]\".*nl.*?\".*)COLORPLACEHOLDER/$${1}#FFCCFF/' \
		| perl -p -e 's/(.* : [[]\".*fr.*?\".*)COLORPLACEHOLDER/$${1}#FFBBDD/' \
		| perl -p -e 's/(.* : [[]\".*fr.*?\".*)COLORPLACEHOLDER/$${1}#FFDDDD/' \
		| perl -p -e 's/(.* : [[]\".*es.*?\".*)COLORPLACEHOLDER/$${1}#FFFFCC/' \
		| perl -p -e 's/(.* : [[]\".*grc.*?\".*)COLORPLACEHOLDER/$${1}#DDFFDD/' \
		| perl -p -e 's/(.* : [[]\".*la.*?\".*)COLORPLACEHOLDER/$${1}#FFFFCC/' \
		| perl -p -e 's/COLORPLACEHOLDER/#CCCCCC/'	\
			> $@
