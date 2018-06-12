#!/usr/bin/env perl

use strict;
use warnings;

use File::stat;

my $startDirectory = $ARGV[0];

my $files = `ls -A -1`;

my @farr = split("\n", $files);

my %gids = ();
my %uids = ();

#my %uidToName = ();
#my %gidToName = ();

my %userTotals = (); # a hash of total file size and number of files by user
my %groupTotals = ();

my $numUnreadable = 0;

my $bytesToGigDivide = (1024 * 1024 * 1024);
dirTraverse($startDirectory);

printStats();

sub printStats {
    my $OUT = *STDOUT;
    
    foreach my $uid (keys(%uids)) {
	my $bytes = $userTotals{$uid}{size};
	my $gb = sprintf("%.2f", $bytes/$bytesToGigDivide);
	
	print $OUT "User $uid:\n";
	print $OUT "\t$gb GB ($bytes bytes)\n";
	print $OUT "\t$userTotals{$uid}{numFiles} files (does not include folders/links)\n";
    }

    foreach my $gid (keys(%gids)) {
	my $bytes = $groupTotals{$gid}{size};
	my $gb = sprintf("%.2f", $bytes/$bytesToGigDivide);
	
	print $OUT "Group $gid:\n";
	print $OUT "\t$gb GB ($bytes bytes)\n";
	print $OUT "\t$groupTotals{$gid}{numFiles} files (does not include folders/links)\n";
    }

    print $OUT "$numUnreadable directories were unreadable.\n";
}


sub dirTraverse {
    my ($theDirectory) = @_;
    my @filesHere = split("\n", `ls -1 '$theDirectory'`);

    foreach my $file (@filesHere) {
	
	$file = $theDirectory . '/' . $file;
	
	if ($file eq '.' or $file eq '..') { next; }
	
	my $stats = stat($file);
	if (!defined($stats)) {
	    print STDERR "Couldn't get the info for $file...\n";
	    next;
	}
	
	my $uid = ${stats}->uid;
	my $uName = `ls -ld '$file' | awk '{print \$3}'`;
	chomp($uName);
#	if (defined($uidToName{$uid})) {
#	    if ($uidToName{$uid} ne $uName) {
#		print STDERR "WARNING: The user ID $uid used to be $uidToName{$uid}, but we just found it used for user $uName as well.\n";
#	    }
#	} else {
#	    $uidToName{$uid} = $uName;
#	}

	my $gid = ${stats}->gid;
	my $gName = `ls -ld '$file' | awk '{print \$4}'`;
	chomp($gName);
#	if (defined($gidToName{$gid})) {
#	    if ($gidToName{$gid} ne $gName) {
#		print STDERR "WARNING: The group ID $gid used to be $gidToName{$gid}, but we just found it used for group $gName as well.\n";
#	    }
#	} else {
#	    $gidToName{$gid} = $gName;
#	}

	if (!defined($userTotals{$uName})) {
	    $userTotals{$uName} = {
		size => 0,
		numFiles => 0,
		};
	}
	if (!defined($groupTotals{$gName})) {
	    $groupTotals{$gName} = {
		size => 0,
		numFiles => 0,
	    };
	}
	
	$uids{$uName} = 1;
	$gids{$gName} = 1;
	
	if (-f $file && !(-l $file)) { # file and not a symlink...
	    # We saw a file, so increment the file count...
	    $userTotals{$uName}{numFiles}++;
	    $groupTotals{$gName}{numFiles}++;
	}
	
	$userTotals{$uName}{size}  += ${stats}->size;
	$groupTotals{$gName}{size} += ${stats}->size;
	
	#print STDERR "$file\n";

	if (-d $file && !(-l $file)) { # file is a directory and NOT a symlink
	    if (!(-x $file) || !(-r $file)) { #  Can't execute the directory, so we can't read it
		$numUnreadable++;
		#print STDERR "Can't read directory $file...\n";
	    }
	    print STDERR "Traversing $file...\n";
	    dirTraverse($file);
	} else {
	    # print "Not traversing $file...\n";
	}
    }
}
