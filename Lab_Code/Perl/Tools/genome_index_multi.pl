#!/usr/bin/perl -w

#@COMMENT@ genom_index_multi.pl will do the indexing to make STAR, Tophat/Bowtie, and other genome indices from an input .fa and .gtf file (both are required). Frequency-of-use rating: 1/10.

# By Alex Williams, 2017

use strict; use warnings; use Getopt::Long;
#use File::Slurp; # <-- for reading an entire file into memory.
# If the system doesn't have slurp, then:
#  sudo perl -MCPAN -e shell
#  install File::Slurp

use Term::ANSIColor;
$| = 1;  # Flush output to STDOUT immediately.

my $isDebugging = 0;
my $verbose = 1; # use 'q' (quiet) to suppress this

sub quitWithUsageError($) { print STDOUT ($_[0] . "\n"); printUsageAndContinue(); warn($_[0] . "\n"); exit(1); }
sub printUsageAndQuit() { printUsageAndContinue(); exit(1); }
sub printUsageAndContinue() {    print STDOUT <DATA>; }
sub debugPrint($) {   ($isDebugging) and print STDERR $_[0]; }
sub verboseWarnPrint($) { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "yellow on_blue"); } }
sub verboseUpdatePrint($) { if ($verbose) { print STDERR Term::ANSIColor::colored($_[0] . "\n", "black on_green"); } }
sub argDie($) { printUsageAndContinue(); die "\n[ERROR: problem in the command line arguments]: $_[0]\n"; }

sub systemZero($) {
	my ($cmd) = @_;
	chomp($cmd); # remove superfluous newline
	print("# Running this command: \n   $cmd\n");
	my $exitCode = system($cmd);
	(0 == $exitCode) or die "ERROR: got the exit code '$exitCode' from running the command above. This is a failure we cannot recover from! The command was: $cmd \n ";
}

my ($gtf, $fasta, $name);
my $santaCruzifyChrNames = 0; # should we add 'chr' to numeric chromosome names? Turns 

my %shouldMake = ();

my $ncpus = 2; # Default = 2

$Getopt::Long::passthrough = 1; # ignore arguments we don't recognize in GetOptions, and put them in @ARGV
GetOptions("help|?|man" => sub { printUsageAndQuit(); }
	   , "q" => sub { $verbose = 0; }
	   , "n|name=s" => \$name
	   , "g|gtf=s" => \$gtf
	   , "f|fa|fasta=s" => \$fasta
	   , "all!"  => sub { $shouldMake{STAR}=$shouldMake{BWA}=$shouldMake{BOWTIE2}=$shouldMake{PERCHR}=$shouldMake{FAI}=$shouldMake{CELLRANGER}=1; }
	   , "star!" => sub { $shouldMake{STAR}=1; }
	   , "bwa!"  => sub { $shouldMake{BWA}=1; }
	   , "cellranger!"  => sub { $shouldMake{CELLRANGER}=1; }
	   , "tophat|bowtie2!" => sub { $shouldMake{BOWTIE2}=1; }
	   , "bowtie!" => sub { die "We cannot currently make the '.ebwt'-format bowtie1 indexes. Specifiy bowtie2 or tophat if you want the '.bt2' format indexes"; }
	   , "perchr|!" => sub { $shouldMake{PERCHR}=1; die "We cannot currently make the per-chromosome fasta directories like 'chr1.fa, chr2.fa, chrN.fa, etc...'. If for some reason you need those, you will need to DIY it."}
	   , "santacruz|cruz|ucsc!" => \$santaCruzifyChrNames
	   , "threads=i" => \$ncpus
	   # "refflat"
	   , "debug!" => \$isDebugging
    ) or printUsageAndQuit();

my $numUnprocessedArgs = scalar(@ARGV);
($numUnprocessedArgs == 0) or quitWithUsageError("[Error] in arguments! You seem to have some extra text / unrecognized parameters in your command line call (or something else is messed up). You must have NO additional command-line arguments to the program! The unrecognized arguments (which may be misspelled, perhaps?) are the following: " . join(" ", @ARGV) . "\n");

my $numStepsToDo = 0;
for my $k (keys(%shouldMake)) { $numStepsToDo++; }

print $santaCruzifyChrNames . "\n";
if ($santaCruzifyChrNames) { ($numStepsToDo == 0) or argDie(qq{if you want to "Santa Cruz-ify" the chromosome names (i.e., turn '1' into 'chr1', 'MT' into 'chrM', etc... for compatibility with the UCSC Genome Browser), you cannot ALSO specify any other steps! You must ONLY Santa-cruz-ify the GTF and FASTA files."}); }
else {                       ($numStepsToDo >= 1) or argDie(qq{missing a required operation--you must specify AT LEAST ONE operation to perform--for example, --star or --tophat or --all !\n (This is usually solved by adding '--all' to the end of your command.)\n}); }

my $memgb_for_cellranger = 8;

my $starDir = "STAR_index_${name}";

my $SANTA_CRUZ_SUFFIX     = ".chr"; # we add '.chr' to indicate that these files now have a CHROMOSOME name appended
my $BWA_FILE_PREFIX       = "bwa.";
my $CELLRANGER_PREFIX     = "cellranger_reference.";
my $CELLRANGER_TEMPDIR_NAME = "CELLRANGER_TMP_".int(rand()*1e10); # random number
my $outp                  = "./";

(defined($name) && length($name) >= 1) or argDie "missing the required 'name' input!";
(defined($fasta)) or argDie "the input 'fasta' file must be specified!";
(defined($gtf))   or argDie "the input 'gtf' annotation file must be specified!";
($fasta !~ /[.](gzip|gz|zip|bzip2|bz2|Z|xz|tar)$/i) or argDie("The fasta file CANNOT be compressed! It must be an uncompressed plain text file.");
($gtf   !~ /[.](gzip|gz|zip|bzip2|bz2|Z|xz|tar)$/i) or argDie("The gtf file CANNOT be compressed! It must be an uncompressed plain text file.");
(-e $fasta) or argDie "the input 'fasta' file must actually exist on the filesystem! We failed to find it at: <$fasta>";
(-e $gtf)   or argDie "the input 'gtf' file must actually exist on the filesystem! We failed to find it at: <$gtf>";

use Cwd 'abs_path';

my $fastaNoSuffix = $fasta; $fastaNoSuffix =~ s/[.](fa|fas|fasta).*$//; # remove any suffix after .fa or .fasta

($name =~ m/[-_.A-Za-z0-9]+/) or argDie "The name must be a SIMPLE name with no weird characters or spaces. Hyphens, dots, and underscores are OK.";
(length($name) < 150)  or argDie "Um... your 'name' for the genome is unreasonably long. Shorten it.";

if ($santaCruzifyChrNames) {
	my $cruzFA  = "${name}.fa";  	#$cruzFA = "${name}.fa #s/[.](fa|fasta)/$SANTA_CRUZ_SUFFIX.$1/i;
	my $cruzGTF = "${name}.gtf"; #gtf}; 	#$cruzGTF     =~ s/[.](gtf)/$SANTA_CRUZ_SUFFIX.$1/i;
	(not -e $cruzFA)  or die "ERROR: File already exists--We are REFUSING to overwrite the file '$cruzFA'!  Exiting early. Delete this file if you want to remake it.";
	(not -e $cruzGTF) or die "ERROR: File already exists--we are REFUSING to overwrite the file '$cruzGTF'! Exiting early. Delete this file if you want to remake it.";
	# for a FASTA file
	my $ucscifyFA  = qq{cat "$fasta" | perl -pe 's/^>\(MT|[MWXYZ]|[0-9]+\)\\b/>chr\$1/' | perl -pe 's/^>chrMT\\b/>chrM/' > $cruzFA};
	systemZero($ucscifyFA);         (-e $cruzFA) or die "Failed to generate the santa-cruz-ified FASTA file.";
	# For a GTF, add "chr" and fix the 'MT' issue
	my $ucscifyGTF = qq{cat "$gtf"   | perl -pe 's/^\(MT|[MWXYZ]|[0-9]+\)\\b/chr\$1/'   | perl -pe 's/^chrMT\\b/chrM/'   > $cruzGTF};
	systemZero($ucscifyGTF);        (-e $cruzGTF) or die "Failed to generate the santa-cruz-ified GTF file.";
	print STDOUT "Looks like we SUCCESSFULLY created the 'Santa Cruz format' output files '$cruzFA' and '$cruzGTF'.\n";
	print STDOUT "Quitting now that this operation has completed...\n[DONE]\n";
	exit(0);
}


my $cellranger_full_path_fasta = abs_path($fasta); # get the ABSOLUTE FULL PATH. This way it will still work if we have to change directories (ugh)
my $cellranger_full_path_gtf   = abs_path($gtf);
(-e $cellranger_full_path_fasta) or die "the input 'fasta' file must actually exist on the filesystem!";
(-e $cellranger_full_path_gtf)   or die "the input 'gtf' file must actually exist on the filesystem!";

my %DAT = ( STAR=>{cmd=>qq{mkdir -p $starDir && STAR --runThreadN $ncpus --runMode genomeGenerate   --genomeDir $starDir   --sjdbGTFfile ${gtf}   --genomeFastaFiles ${fasta}  --sjdbOverhang 100}
		   , outputs=>[qq{${starDir}/genomeParameters.txt}]
		   , reqExes=>["STAR"]
		  }
#STAR_index_$(SPECIES)/genomeParameters.txt: $(SPECIES).fa $(SPECIES).gtf
#	@echo "Creating STAR aligner index output directory <"$(dir $@)">"
#	@echo "This takes about 3 hours on our server for a human / mouse genome. It is substantially faster than a bowtie2 index."
#	@echo "Total file size is about 10x the size of the input FQ file."
#	mkdir -p $(dir $@)
#	STAR  --runThreadN 6   --runMode genomeGenerate   --genomeDir $(dir $@)   --sjdbGTFfile $(SPECIES).gtf   --genomeFastaFiles $(SPECIES).fa   --sjdbOverhang 100
#	@echo "Note from Alex about STAR parameters: From https://www.biostars.org/p/93883/:  --sjdbOverhang is used only at the genome generation step, and tells STAR how many bases to concatenate from donor and acceptor sides of the junctions. If you have 100b reads, the ideal value of --sjdbOverhang is 99, which allows the 100b read to map 99b on one side, 1b on the other side. One can think of --sjdbOverhang as the maximum possible overhang for your reads. (Note: --sjdbOverhang 100 is the default value as of May 2016)"
#	@echo "Note from Alex about STAR parameters: From https://www.biostars.org/p/93883/:  On the other hand, --alignSJDBoverhangMin is used at the mapping step to define the minimum allowed overhang over splice junctions. For example, the default value of 3 would prohibit overhangs of 1b or 2b."
#	@echo "CONCLUSION: We should now have generated the followng eight files in the index subdirectory: 1. Genome  2. SA  3. SAindex  4. chrLength.txt  5. chrName.txt  6. chrNameLength.txt  7. chrStart.txt  8. genomeParameters.txt"

	    , BOWTIE2=>{cmd=>qq{bowtie2-build "$fasta" "${outp}${name}" }
			 , outputs=>["${outp}${name}.1.bt2","${outp}${name}.rev.1.bt2"]
			 , reqExes=>["bowtie2-build"]
			}
	    , CELLRANGER=>{cmd=>(qq{mkdir -p "${CELLRANGER_TEMPDIR_NAME}" && cd "${CELLRANGER_TEMPDIR_NAME}" }
				 . qq{ && }
				 . qq{ cellranger mkref --nthreads=$ncpus --memgb=$memgb_for_cellranger --genome="${name}" --fasta="$cellranger_full_path_fasta" --genes="$cellranger_full_path_gtf" }
				 . qq{ && } . qq{ cd ../ }
				 . qq{ && } . qq{ mv "${CELLRANGER_TEMPDIR_NAME}/${CELLRANGER_PREFIX}${name}" "${outp}${CELLRANGER_PREFIX}${name}" }
				 . qq{ && } . qq{ mv "${CELLRANGER_TEMPDIR_NAME}/Log.out"                     "${outp}${CELLRANGER_PREFIX}${name}/Log.out" }
				 . qq{ && } . qq{ rmdir "${CELLRANGER_TEMPDIR_NAME}" }
				 . qq{ })
			   , outputs=>["${outp}{CELLRANGER_PREFIX}${name}/reference.json"]
			   , reqExes=>["cellranger"]}
	    , BWA=>{cmd=>qq{bwa index -p ${outp}${BWA_FILE_PREFIX}${name}  $fasta }
		    , outputs=>["${outp}${BWA_FILE_PREFIX}$name.amb","${outp}${BWA_FILE_PREFIX}$name.ann","${outp}${BWA_FILE_PREFIX}$name.pac"]
		    , reqExes=>["bwa"]
		   }
	    , FAI=>{cmd=>qq{samtools faidx $fasta && TTT=`mktemp` && mv -f ${fasta}.fai "\$TTT" && mv "\$TTT" ${outp}${fasta}.fai && ln -f -s ${outp}${fasta}.fai ${outp}${fastaNoSuffix}.fai && cut -f 1,2 ${outp}${fasta}.fai | sort -k1,1 > ${outp}chrom.sizes.txt}
		    , outputs=>["${outp}${fasta}.fai", "${outp}${fastaNoSuffix}.fai","${outp}chrom.sizes.txt"]
		    , reqExes=>["samtools", "cut", "sort"]
		   }
	    
## Generates the split-up 'chr1.fa, chr2.fa... etc...' files instead of just 'genome.fa' in one single file.
#$(SPECIES)_chromosome_fastas/makefile_ran_successfully.touch.txt: INPUT_PER_CHROM_FASTA_DIR
#	mkdir -p $(dir $@)
#	for f in $</*.fa; do \
#		base=`basename $$f` ; \
#		if [[ "$$base" == "MT.fa" ]] ; \
#			then base="M.fa" ; \
#		fi ; \
#		out="$(dir $@)/chr$$base" ; \
#		echo "Generating the proper single-chromosome fasta file for input file $$f. Output will be $$out..."; \
#		head -n 1  $$f | perl -pe 's/^>MT/>M/' | perl -pe 's/^>/>chr/' > $$out ; \
#		tail -n +2 $$f >> $$out ; \
#	done
#	touch $@
#
#	/bin/mv -f ./Log.out $(dir $@)/

	  );


for my $stepName (keys(%shouldMake)) {
	print $stepName . "\n";
	if (!exists($DAT{$stepName})) {
		print "problem...\n";
		next;
	}
	my $cmd = $DAT{$stepName}{cmd};
	my @outs = @{$DAT{$stepName}{outputs}};

	systemZero($cmd);

	for my $requiredOutFile (@outs) {
		(-e $requiredOutFile) or die "[ERROR]: Failure in the '$stepName' step---we failed to generate the expected output file '$requiredOutFile'! Check the path and capitalization.\n";
	}
}

print("[DONE]");


exit(0); # looks like we were successful


################# END MAIN #############################

__DATA__

genome_index_multi.pl --name=STRING --gtf=FILE --fasta=FILE   OPERATION
genome_index_multi.pl     -n STRING    -g FILE      -f FILE   OPERATION

By Alex Williams

Generates the required genome indexes for various program. See below for specifics.

syntax: genome_index_multi.pl --name=MyGenome --gtf=GenomeAnnot.gtf --fasta=Seqs.fasta --all

Two kinds of inputs are required, PARAMETERS and OPERATIONS.

PARAMETERS (3 required):
     You must always supply --name, --gtf (file), and --fasta (file).
     --name  / -n : Name for output files, like "hippopotamus" or "My_Genome". No spaces or weird characters.
     --gtf   / -g : Input GTF file. Must exist. Cannot be compressed.
     --fasta / -f : Input FASTA file. Must exist. Cannot be compressed.

OPERATIONS (specify which indexes you want generated, or '--all' for all of them)
  --all: Generate ALL the various index files described above! Recommended option.
Or you can only do a subset to generate indexes just for a particular program:
  --star: Generates STAR indexes
  --cellranger: Generates CELLRANGER indexes, which are STAR indexes plus some metadata.
  --bwa: Generates indexes for the BWA aligner.
  --tophat / --bowtie2: Generates a bowtie2-format set of indexes (the ones ending in .bt2).
  --bowtie (version 1): NOT SUPPORTED right now. 

Examples:

Make ALL indexes (for BWA, Bowtie, STAR, Cellranger), name the final genome "Hippopotamus":
    genome_index_multi.pl --name=Hippopotamus --gtf=HippoGenes.gtf --fasta=Hipp2017.fasta --all

Make ONLY the 'STAR' aliger indexes:
    genome_index_multi.pl -n HippoSTAR  -g HGenes.gtf  -f Hip10.fa  --star

Convert an ENSEMBL fasta and gtf file to UCSC format chromosome names (1 -> chr1, etc...):
    genome_index_multi.pl --name=Hippopotamus --gtf=HippoGenes.gtf --fasta=Hipp2017.fasta --santacruz


Other options:

  --perchr: NOT SUPPORTED YET. Would generate individual files for each chromosome.

  --threads: Specify number of processors/CPUs to use. Default = 2.

  --debug: Unused

--santacruz: Reformat the input FASTA and GTF files to have 'UCSC Genome Browser'-compatible names.
  Turns '1' -> 'chr1', '2' -> 'chr2', ... 'MT' -> 'chrM'. It just adds 'chr' to any only-numeric chromosome
  or any chromosome named:  W, X, Y, Z, M or MT. Note that MT is the only two-letter chromosome, and it does not
  fit the pattern of the others, because it becomes 'chrM' (mitochondrial) instead of 'chrMT' (which is not recognized)
  This option must be run ALONE, it does not work with any other generation command.
         index_genomes.pl --name=NewUCSC --fa=A.FA --gtf=A.GTF --santacruz
  Will generate exactly two new files:
         NewUCSC.FA and NewUCSC.GTF
  Requires that the new files do not already exist--will not overwrite existing files.

Bugs / To-do:

None known at the moment.

================
