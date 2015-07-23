#!/usr/bin/perl -w
#use Spreadsheet::ParseExcel;
#use Getopt::Long;
#use strict;

print("\n");
print("*"x80 . "\n");
print("*"x80 . "\n");
print("*"x80 . "\n");
print("excel2tab.pl is deprecated: you should use the following instead:\n");
print("\n");
print("Either:   git clone https://github.com/dilshod/xlsx2csv.git\n");
print(" (should be installable with 'sudo pip install xlsx2scv)\n");
print(" Then you can just type: xlsx2csv -d tab -p 1 INPUT.xlsx OUT.txt\n");
print(" (That sample command only converts the FIRST sheet--omit '-p 1' to get all sheets.)\n");
print("Note: any cell with a comma in it will become QUOTED. Quotes turn into 2x or 3x quotes.\n");
print("\n");
print("Alternate option which requires LibreOffice: unoconv (installable via yum / apt-get).\n");
print("*"x80 . "\n");
print("*"x80 . "\n");
print("*"x80 . "\n");
print("\n");

exit(1); # exit with an error to indicate that this is non-functional
  
# my $directory_suffix = ".sheets";
# my $worksheet_suffix = ".tab";

# my $force_overwrite;
# my $list_mode;
# my $directory_prefix;
# my $print_help;

# GetOptions("f"   => \$force_overwrite,
# 	   "l"   => \$list_mode,
# 	   "h"   => \$print_help);

# if ($print_help) {
#     print STDOUT <DATA>;
#     exit(0);
# }

# sub makeSheetsDirectory {
#     my ($filename) = @_;
#     my $dirname = $filename . $directory_suffix;
#     if (-e $dirname) {warn "$dirname already exists";}
#     if (!(-d $dirname) && !mkdir($dirname)) {
# 	warn "Can't create directory $dirname, skipping $filename";
# 	return 0;
#     }
#     return ($dirname . '/');
# }

# ## The skeleton was taken from the module documentation
# foreach my $file (@ARGV) {
#     my $excel = Spreadsheet::ParseExcel::Workbook->Parse($file);
#     print "File: $file\n" if (scalar(@ARGV) > 1);
#     my $dirprefix;
#     $dirprefix = makeSheetsDirectory($file) unless ($list_mode);
#     next if ($list_mode && !$dirprefix);
#     foreach my $sheet (@{$excel->{Worksheet}}) {
# 	printf("Sheet: %s\n", $sheet->{Name});
# 	next if $list_mode;
# 	my $tabfilename = $dirprefix . $sheet->{Name} . $worksheet_suffix;
# 	open OUT, ">", $tabfilename or die "couldn't open file $tabfilename";
# 	$sheet->{MaxRow} ||= $sheet->{MinRow};
# 	foreach my $row ($sheet->{MinRow} .. $sheet->{MaxRow}) {
# 	    $sheet->{MaxCol} ||= $sheet->{MinCol};
# 	    foreach my $col ($sheet->{MinCol} ..  $sheet->{MaxCol}) {
# 		print OUT "\t" if ($col > $sheet->{MinCol});
# 		my $cell = $sheet->{Cells}[$row][$col];
# 		print OUT $cell->{Val} if ($cell);
# 	    }
# 	    print OUT "\n";
# 	}
# 	close OUT;
#     }
# }

# __DATA__
# excel2tab.pl - parse an Excel .xls file, and print each worksheet in the
# file to a tab-delimited version.  The worksheets from each .xls file will
# be placed in a directory named by appending .sheets to the input file name.
# Each tab file will have a name matching the worksheet it was created from, 
# with a .tab appended.

# Usage:
#     excel2tab.pl [options] file1.xls file2.xls

# Options:

#     -h       print this help message
#     -f       force overwrite of existing files and directories
#     -l       just list the worksheets in each file, do not write output
