#!/usr/bin/perl -w


# trash.pl is meant to be a "safer rm," which moves things to a trash can
# instead of immediately removing them. trash.pl will mangle names, so it's
# not necessarily going to be easy to restore things.

# Any slashes will be replaced with __PATH__, which means that restoring
# any directory structure will be exceedingly annoying. So be careful!

#
# Use and abuse as much as you want.
# Put it in /usr/bin/ or $HOME/bin/
# Original by Daniel Cote
#
# Modified by Alex Williams
#
# Most recent version of this file available at
# http://www.novajo.ca/trash.html
#
# Instead of deleting files, this script moves them to the appropriate trash folder.
# Trash folders are either private ($HOME/.Trash/) or public (/Volumes/name/.Trashes/
# and /.Trashes).  In the directory .Trashes/, each user has a trash folder named after
# its userid (501, 502, etc...).
#
# This script will simply stop if anything goes wrong. If by any chance it appears
# that it might overwrite a file, it will let you know and ask for a confirmation.
#
# This script could be used as a replacement for rm.  The options passed as arguments
# are simply discarded.

# Usage: trash [options] file1|dir1 file2|dir2 file3|dir3 ...
# You can use wildcards. To erase a directory, you are probably better off just naming
# the directory (e.g. trash Documents/) as opposed to trash Documents/* since the first
# case will actually keep the hierarchy and simply move Documents into the trash whereas
# the second one will move each item of the document folder into the trash individually.
#

use strict;
use warnings;

use Cwd ('abs_path', 'getcwd');
use File::Basename;

use Term::ANSIColor;

my $outputIsColorTerminal = (-t STDOUT);

# sub hasBashColor() {
#     {
# 	## Check that the TERM is a color one...
# 	my $theTerm = $ENV{'TERM'};
# 	if (!defined($theTerm)) { return 0; }
# 	if (!(($theTerm eq 'xterm-color') || ($theTerm eq 'xterm-256color'))) { return 0; }
#     }

#     ## Check that ths shell is BASH, because nothing else will understand these escape sequences...
#     my $theShell = $ENV{'SHELL'};
#     if (!defined($theShell)) { return 0; }
    
#     if ($theShell =~ /.*bash$/) { return 1; }
#     else { return 0; }

# }

sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

if ( scalar @ARGV == 0) {
    die "trash.pl needs at least one argument. It works in a similar fashion to *rm*.\nUsage: trash.pl thing1 thing2 thing3 ...\n\n";
}

my $MAX_INDEX = 1000; # <-- number of duplicate-named files allowed in a directory

my $username = `whoami`;
chomp($username);
if (!(length($username) > 0)) {
    die qq{We could not get your username! We need it to make the trash directory in /tmp!\n};
}
if (!defined($username) || (length($username) <= 0) || $username =~ /[\/\\"' 	*+]/) {
    die qq{The username is either blank, or it has either a space in it, or a tab, or a star, or a plus, or a quotation mark, or a slash, or a backslash! We cannot reliably create a trash directory without it being dangerous!!!\n};
}





my $tab = " " x length("trash.pl: "); # "tab" is just the number of spaces in the program name

print STDOUT (qq{-} x 80) . "\n";


sub fatalComplaint($$) {
    anyComplaint(@_, "red");
    exit(1);
}

sub minorComplaint($$) {
    return(anyComplaint(@_, "yellow"));
}

sub anyComplaint($$$) {
    my ($aboutWhat, $why, $color) = @_;
    if ($outputIsColorTerminal) { print STDOUT color($color); }
    print ("trash.pl: Not deleting \"$aboutWhat\": " . $why . "\n");
    if ($outputIsColorTerminal) { print STDOUT color("reset"); }
}


foreach my $itemToCheck (@ARGV) {
    ## Go through and CHECK first to see what files are super duper unsafe for deletion, and quit if we find any.
    $itemToCheck = trim($itemToCheck);

    my $HOME = $ENV{"HOME"};
    if ($itemToCheck =~ /^[\/]*$HOME[.\/]*$/) {
	fatalComplaint($itemToCheck, "you cannot remove the home directory! You can use /bin/rm to attempt to remove it if you really must.")
    }

    if (($itemToCheck =~ /^[.\/]*Applications[.\/]+Utilities[.\/]*$/)
	or ($itemToCheck =~ /^[.\/]*(System|Developer).*/)  ## protect *any* subdirectory of these
	or ($itemToCheck =~ /^[.\/]*$HOME[.\/]+System.*/) ## protect *any* subdirectory of these
	or ($itemToCheck =~ /^[.\/]*$HOME[.\/]+(Desktop|Library|Documents)[.\/]*$/) ## protect only the ROOT level of these
	or ($itemToCheck =~ /^[.\/]*(Applications|sbin|bin|boot|dev|etc|lib|mnt|opt|sys|usr|var)[.\/]*$/) ## protect only the ROOT level of these
	) {
	fatalComplaint($itemToCheck, "removing this system file is unsafe! Trash.pl will not peform this operation! You can use /bin/rm if you REALLY want to remove it. Quitting.");
    }
    
    if (($itemToCheck =~ /^[.\/]+$/) ## any combination of just dots or slashes, and nothing else
	) {
	fatalComplaint($itemToCheck, "removing a file that is the current directory or is directly above the current directory is not something that trash.pl thinks is safe! Trash.pl will not peform this operation! You can use /bin/rm if you REALLY want to remove it. Quitting.");
    }
}

my $numItemsDeleted = 0;
foreach my $itemToDelete (@ARGV) {
    #print $itemToDelete . "\n";
    $itemToDelete = trim($itemToDelete);
    
    if ($itemToDelete =~ m|^-| && (!(-l $itemToDelete)) && (!(-e $itemToDelete))) {
	# If the option starts with a hyphen AND it isn't a file or symlink
	minorComplaint($itemToDelete, "it was not a filename, and it appears to have been intended as a command line switch, which trash.pl does not make any use of.");
    }
    
    if ((not (-l $itemToDelete)) and (not (-e $itemToDelete))) {
	# The item doesn't even exist (and it isn't a symlink), so skip it
	minorComplaint($itemToDelete, "it did either not exist or could not be read. Check that this file actually exists.");
	next;
    }
	
    if (not (-f $itemToDelete || -d $itemToDelete || -l $itemToDelete)) {
	# The thing to delete has to be either a file (-f), a directory (-d), or a symlink (-l).
	minorComplaint($itemToDelete, "it was not a file, directory, or symlink. You can use /bin/rm to remove it if you really want to.");
	next;
    }

    my $thepath = undef;

    if (-l $itemToDelete) {
	# the "item to delete" is a symbolic link. Therefore, do NOT follow the symlink and delete the real file. Instead, just delete the symlink.
	$thepath = $itemToDelete; # <-- just delete the symlink, DON'T follow the link to the real file!!!
    } else {
	$thepath = abs_path($itemToDelete); # <-- figure out what the full path of this file is.
    }

    #print "Item to delete: " . $thepath . "\n";

    my $basename = basename($thepath);
    my $dirname  = dirname($thepath);

    #print "Full name of the thing we are deleting to make: " . $thepath . "\n";
    #print "Location to put it: " . "/tmp/alexgw/Trash" . $thepath . "\n";

    #print "Specific file/folder: " . $basename . "\n";
    #print "Directory that it is in: " . $dirname . "\n";
    #print "\n";

    my $trash = qq{/tmp/${username}/Trash};

    if (! -e $trash) {
	print qq{trash.pl: Making trash directory "$trash"\n};
	`mkdir -p $trash`;
	`chmod og-rwx $trash`; # making it so no one else can read it...
	`chmod u+rwx $trash`; # but we can read it...
    }

    if (!(-d $trash)) {
	die "Error: $trash is not a directory";
    }

    if (!(-l $thepath) && !(-e $thepath)) {
	# need to check that there is a symlink (-l) at the location, OR a real thing (-e)
	die "Error getting full path to file: $thepath does not exist\n";
    }

    my $UNSAFE_CHARS = q{"'*+\$;};
    if ($thepath =~ /[$UNSAFE_CHARS]/) {
	if ($outputIsColorTerminal) { print STDERR color("red"); }
	print STDERR "trash.pl: Uh oh, trash.pl has no idea how to deal with filenames with any of these " . length($UNSAFE_CHARS) . " unusual characters in them: $UNSAFE_CHARS . That might be bad news, so we are just going to abort. Try using the real version of 'rm' in /bin/rm in this case.\n";
	if ($outputIsColorTerminal) { print STDERR color("reset"); }
	exit(1)
    }

    my $index = 1;
    my $trashDupeSuffix = qq{.trash_duplicate};

    my $trashedDirLoc  = qq{$trash/$dirname};

    if (! -d $trashedDirLoc) {
	#print qq{trash.pl: Making a new trash subdirectory at $trashedDirLoc\n};
	system(qq{mkdir -p "$trashedDirLoc"});
    }

    my $trashedFileLoc =
	($thepath =~ m{^\/})
	? "${trash}${thepath}"
	: "${trash}/${thepath}";

    if (-e $trashedFileLoc) {
	my $firstFailedAttempt = $trashedFileLoc;
	while (-e $trashedFileLoc) {
	    $index++;
	    $trashedFileLoc =
		($thepath =~ m{^\/})
		? "${trash}${thepath}${trashDupeSuffix}.${index}"
		: "${trash}/${thepath}${trashDupeSuffix}.${index}";
	    
	    if ($index >= $MAX_INDEX) {
		if ($outputIsColorTerminal) { print STDERR color("red"); }
		print STDERR qq{Error trying to rename \"$itemToDelete\" to avoid overwriting an existing file with the same name that is already in the trash (there are too many files with the same name).\nYou should probably empty the trash.};
		if ($outputIsColorTerminal) { print STDERR color("reset"); }
		exit(1)
	    }
	}

	if ($outputIsColorTerminal) { print STDOUT color("yellow"); }
	print STDOUT "trash.pl: There was already a file in the trash named\n";
	print STDOUT qq{${tab}${firstFailedAttempt}\n};
	print STDOUT qq{${tab}so \"};
	print STDOUT qq{${trashDupeSuffix}.${index}};
	print STDOUT qq{\" was appended to this filename before it was trashed.\n};
	if ($outputIsColorTerminal) { print STDOUT color("reset"); }
    }

    my $mvResult = system(qq{mv "$thepath" "$trashedFileLoc"});
    
    if (0 != $mvResult) {
	## mv operation FAILED for some reason...
	if ($outputIsColorTerminal) { print STDERR color("red"); }
	print STDERR "trash.pl: Failed to trash the file $itemToDelete.\n          Probably you do not have the permissions required to move this file.\n";
	if ($outputIsColorTerminal) { print STDERR color("reset"); }
    } else {
	## Move succeeded...
	if (-l $itemToDelete) {
	    if ($outputIsColorTerminal) { print STDOUT color("green"); }
	    print STDOUT qq{trash.pl: Trashing the symbolic link $itemToDelete -> $trashedFileLoc\n};
	    if ($outputIsColorTerminal) { print STDOUT color("reset"); }
	} else {
	    if ($outputIsColorTerminal) { print STDOUT color("green"); }
	    print STDOUT qq{trash.pl: };
	    print STDOUT qq{$itemToDelete};
	    print STDOUT qq{ -> };
	    print STDOUT qq{$trashedFileLoc\n};
	    if ($outputIsColorTerminal) { print STDOUT color("reset"); }
	}
    }
    $numItemsDeleted++;
}

if ($numItemsDeleted > 0) {

    print STDOUT qq{trash.pl: Note: the trash directory will be auto-deleted every reboot, without warning.\n};

    print STDOUT (qq{-} x 80) . "\n";

}
