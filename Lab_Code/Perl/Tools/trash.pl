#!/usr/bin/perl -w

#@COMMENT@ trash.pl is a "safer rm" that moves files to a trash directory in /tmp/ instead of immediately removing them. trash.pl will mangle names, so it's not necessarily going to be easy to restore things. But it is probably better than instantly deleting them. May fill up your /tmp partition if you delete extremely large files. Any slashes in a filename will be replaced with __PATH__, which means that restoring any directory structure will be exceedingly annoying. So be careful! Frequency-of-use rating: 1/10

# Use and abuse as much as you want.
# Put it in /usr/bin/ or $HOME/bin/
# Original code by Daniel Cote
# Modified by Alex Williams, 2015
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
#use threads;

my $outputIsColorTerminal = (-t STDOUT); # Check if the output is to the TERMINAL (then we should print color), or to a redirect (then don't print color!)

my $debugMode = 0;
my $pn = "trash.pl"; # program name

use File::Path qw(make_path); # make_path is basically "mkdir -p"
use Cwd ('abs_path', 'getcwd');
use File::Copy; # copy(...) and move(...)

eval {
    require File::Basename;
    File::Basename->import();
}; # <-- mandatory semicolon!!! or it breaks horribly
if($@) {
    #$GLOBAL_HAS_BASENAME = 0;
    warn "Note: File::Basename Perl library is not installed / not in the include path.";
    exit(1);
    # FAILURE!!! Cannot find File::Basename apparently
}

eval {
    require Term::ANSIColor;
    Term::ANSIColor->import();
}; # <-- mandatory semicolon!!! or it breaks horribly
if($@) {
    $outputIsColorTerminal = 0;
    warn "Note: Term::ANSIColor Perl library is not installed / not in the include path. Disabling color output.";
    # FAILURE!!! Cannot find Term::ANSIColor apparently
}

sub trimWhitespace($) { my $str = shift; $str =~ s/^\s+//; $str =~ s/\s+$//; return $str; }
sub getFilesystemIDForPath($) {
    my $path = shift;
    if (-l $path) { return ((lstat($path))[0]); } # lstat is for SYMLINKS
    else {          return ((stat($path))[0]);  } # "stat" is a perl built-in. Returns a number! Doesn't work on broken symlinks.
}

my $MAX_DUPE_FILENAMES = 1000; # <-- number of duplicate-named files allowed in a directory before we just give up

my $USER = getlogin() || ((getpwuid($<))[0]) || undef; # Note the fallbacks---getlogin SHOULD work however!
(defined($USER) && (length($USER)>0)) or die qq{We could not get your current logged-in username! We need it to make the trash directory in /tmp!\n};
($USER =~ m/^[-_.\w]+$/)              or die qq{The username (\"$USER\") contains DANGEROUS characters in it that would not be suitable for making a directory! (Examples: spaces, or a star, or a plus, or a quotation mark, or a slash, or a backslash, etc etc! We cannot reliably create a trash directory without it being dangerous!!! Quitting.\n};

my $TRASH_ROOT = qq{/tmp/$USER/Trash}; # <---------- HARD COED

my $tab = " " x length("$pn: "); # "tab" is just the number of spaces in the program name

sub printComplaintAboutFile($$$) {
    my ($aboutWhat, $why, $color) = @_;
    chomp($why);
    if ($outputIsColorTerminal) { print STDOUT color($color); }
    print STDOUT ("$pn: Not deleting \"$aboutWhat\": " . $why . "\n");
    if ($outputIsColorTerminal) { print STDOUT color("reset"); }
}

sub fatalExit($$)   { printComplaintAboutFile($_[0], $_[1], "red"); exit(1); }
sub complainAboutFile($$) { return(printComplaintAboutFile($_[0], $_[1], "yellow")); }

sub printInfo($) {
    my ($info) = @_;
    chomp($info);
    if ($outputIsColorTerminal) { print STDOUT color("green"); }
    print STDOUT ("$pn: $info\n");
    if ($outputIsColorTerminal) { print STDOUT color("reset"); }
}

# ========================================================================

(scalar(@ARGV) > 0) or die "$pn needs at least one argument--the file(s) to delete. It works in a similar fashion to *rm*.\nUsage: $pn file1 file2 file3 ...\n\n";

print STDOUT (qq{-} x 80) . "\n";

if (not (-d ${TRASH_ROOT})) {
    printInfo(qq{Making new trash directory: \"${TRASH_ROOT}\"});
    File::Path::make_path(${TRASH_ROOT}, { verbose => ($debugMode ? 1 : 0), mode => 0700 } ) or die "Failed to make trash root directory \"${TRASH_ROOT}\". This may be a permissions issue."; # 0700 = we can read/write this file, but no one ELSE can!
    (-d ${TRASH_ROOT}) or die "Failed to create trash root directory.";
}

my $TRASH_ROOT_FILESYSTEM_ID = getFilesystemIDForPath($TRASH_ROOT); # The actual physical(?) filesystem ID. We will use this check to see if a "mv" operation just moves the file on the same filesystem (fast) or has to copy it to the trash on a DIFFERENT filesystem (slow!)

foreach my $itemToCheck (@ARGV) {
    ## Go through and CHECK first to see what files are super duper unsafe for deletion, and quit if we find any.
    $itemToCheck = trimWhitespace($itemToCheck);
    my $HOME = $ENV{"HOME"};
    if ($itemToCheck =~ /^[\/]*$HOME[.\/]*$/) {
	fatalExit($itemToCheck, "you cannot remove the home directory! You can use /bin/rm to attempt to remove it if you really must.")
    }
    if (($itemToCheck =~ /^[.\/]*Applications[.\/]+Utilities[.\/]*$/i)
	or ($itemToCheck =~ /^[.\/]*(System|Developer).*/i)  ## protect *any* subdirectory of these
	or ($itemToCheck =~ /^[.\/]*$HOME[.\/]+System.*/i) ## protect *any* subdirectory of these
	or ($itemToCheck =~ /^[.\/]*$HOME[.\/]+(Desktop|Library|Documents)[.\/]*$/i)
	or ($itemToCheck =~ /^[.\/]*(Applications|sbin|bin|boot|dev|etc|lib|mnt|opt|sys|usr|var)[.\/]*$/i)
	or ($itemToCheck =~ /^[-]*(fr|rf|r|f)$/i) ## assume that r / f in either combination, with our without the hyphen(s) is a mistake since trash.pl doesn't take -rf operators.
	) {
	fatalExit($itemToCheck, "removing this file seems potentially unsafe! $pn will not peform this operation! You can use /bin/rm if you REALLY want to remove it. Quitting.");
    }    
    if (($itemToCheck =~ /^[.\/\\]+$/)) { ## any combination of just dots or slashes/backslashes, and nothing else
	fatalExit($itemToCheck, "removing a file that is the current directory or is directly above the current directory is not something that $pn thinks is safe! $pn will not peform this operation! You can use /bin/rm if you REALLY want to remove it. Quitting.");
    }
}

my $numItemsDeleted = 0;
foreach my $darg (@ARGV) {
    $darg = trimWhitespace($darg); # darg is the argument to delete
    if (($darg =~ m|/+$|) and (not -l $darg) and (not -d $darg)) {
	fatalExit($darg, "ends in a trailing slash, but was not a symlink or a directory!");
	next;
    }
    $darg =~ s|/+$||; # Remove any extraneous trailing slash(es) on any item to delete. Solves the problem of us deleting directories by following a symlink! Important; otherwise "rm SYMLINK/" removes the real directory instead of the symlink!
    my $isSymlink = (-l $darg);

    if ($darg =~ m/^-/ && !$isSymlink && !(-e $darg)) { # If the argument starts with a hyphen AND it isn't a file or symlink, it's probably a mistaken attempt at a command line argument.
	complainAboutFile($darg, "it was not a filename, and it appears to have been intended as a command line switch, which $pn does not make any use of.");
	next;
    }
    
    if (!$isSymlink and !(-e $darg)) {
	complainAboutFile($darg, "it did either not exist (and is not a symlink) or could not be read. Check that this file actually exists.");
	next;
    }
	
    if (not (-f $darg or -d $darg or $isSymlink)) { # The thing to delete has to be either a file (-f), a directory (-d), or a symlink (-l).
	complainAboutFile($darg, "it was not a file, directory, or symlink. You can use /bin/rm to remove it if you really want to.");
	next;
    }

    my $delPath = $isSymlink ? $darg : abs_path($darg);   # <-- If it's a symlink (-l), then DO NOT FOLLOW the symlink to the real file---just delete the symlink! Otherwise get the abs_path
    my $basename = File::Basename::basename($delPath);
    my $dirname  = File::Basename::dirname($delPath);

    #print "Full name of the thing we are deleting to make: " . $delPath . "\n";
    #print "Location to put it: " . "/tmp/alexgw/Trash" . $delPath . "\n";

    #print "Specific file/folder: " . $basename . "\n";
    #print "Directory that it is in: " . $dirname . "\n";
    #print "\n";
    
    my $fsys = getFilesystemIDForPath($darg);
    defined($fsys) or fatalExit($delPath, "$fsys was not defined, probably this is some issue with it being a file type we didn't anticipate. Maybe it's a symlink and my code for handling symlinks is bad.\n");
    my $moveIsOnSameFilesystem = ($fsys == $TRASH_ROOT_FILESYSTEM_ID); # Is our 'mv' operation actually COPYING the file? (To a different /tmp filesystem?)

    ($isSymlink or ((-e $delPath) and (-r $delPath))) or fatalExit($delPath, "it was not possible to get the full path to this file, as it appears not to exist (or possibly is unreadable). For non-symlinks, we REQUIRE that the file exists and is readable if we are going to try to throw it away.");

    #                vvvv This list below is of the SAFE characters. \w is alphanumeric "word" characters.
    ($delPath =~ m/^[-\/\\\w .,;:\!\@\#\&\(\)\{\}\[\]=]+$/) or fatalExit($delPath, "it had certain 'unsafe' characters that we did not know how to handle. That might be bad news, so we are just going to abort. Try using the real version of 'rm' in /bin/rm in this case.\n");

    my $index = 1;
    my $trashDupeSuffix = qq{.trash_dupe};
    my $trashedDirLoc   = qq{${TRASH_ROOT}/$dirname};

    (-d $trashedDirLoc) or File::Path::make_path($trashedDirLoc, { verbose => ($debugMode ? 1 : 0), mode => 0700 } ) or die "Failed to make \"$trashedDirLoc\".";

    my $moveTo = "${TRASH_ROOT}/${delPath}"; # Always a slash, even if there was one already
    if (-e $moveTo) {
	# If the file ALREADY EXISTS at the target location
	my $failedMoveAttempt = $moveTo;
	while (-e $moveTo) {
	    $index++;
	    $moveTo = "${TRASH_ROOT}/${delPath}${trashDupeSuffix}.${index}";
	    ($index < $MAX_DUPE_FILENAMES) or fatalExit($delPath, "Error trying to rename this file to avoid overwriting an existing file with the same name that is already in the trash (there are too many files with the same name).\nYou should probably empty the trash.");
	}

	if ($outputIsColorTerminal) { print STDOUT color("yellow"); }
	print STDOUT "$pn: There was already a file in the trash named\n";
	print STDOUT qq{${tab}${failedMoveAttempt}\n};
	print STDOUT qq{${tab}so \"};
	print STDOUT qq{${trashDupeSuffix}.${index}};
	print STDOUT qq{\" was appended to this filename before it was trashed.\n};
	if ($outputIsColorTerminal) { print STDOUT color("reset"); }
    }


    my $mvOK; # Was the move operation OK (successful)?
    my $isBackgroundMove = 0;

    if ($moveIsOnSameFilesystem) {
	$mvOK = File::Copy::move($delPath, $moveTo); # Just MOVE the file, nothing fancy.
    } else {
	# Move operation CROSSES FILESYSTEMS, which means it will be SLOW...
	# ...so what we do is that we run a shell background task to "mv file1 trash/ &" (with the '&').
	# This will fail if the user closes the terminal while the deletion is still in progress (which will interrupt the deletion)
	# but is otherwise robust and allows "trash.pl" to exit while the job continues.
	# The file is named old_filename.DELETION_IN_PROGRESS until the delete operation succeeds, at which point it is removed.
	my $delTEMP = $delPath . ".DELETION_IN_PROGRESS";
	(not -e $delTEMP)                    or fatalExit($delPath, "the $delTEMP temp file already exists... this may be a result of a failed 'rm' operation earlier. Try using the system '/bin/rm FILENAME' to delete the file instead.");
	File::Copy::move($delPath, $delTEMP) or fatalExit($delPath, "The command 'move' failed when moving from '$delPath' to '$delTEMP'. Try using the system '/bin/rm FILENAME' to delete the file instead.");
	system(qq{/bin/mv "$delTEMP" "$moveTo" &}); # move in the ******BACKGROUND********. Note that we CANNOT get the exit code from this, because the job runs in the background (it doesn't finish here)
	$mvOK = 1; # We have to just assume it worked, since we can't actually tell if it succeeded (since the move happens in the background, and might finish an hour from now)
	$isBackgroundMove = 1;
    }

    if ($mvOK) {
	my $prefix = $isSymlink        ? qq{Trashing the symlink } : qq{};
	my $suffix = $isBackgroundMove ? qq{ (moving file to a different filesystem)} : q{};
	printInfo(qq{${prefix}${darg} -> ${moveTo}${suffix}});
    } else {
	printComplaintAboutFile($delPath, "Probably you do not have the permissions required to move this file: $!.", "red");
	next;
    }
    $numItemsDeleted++;
}

if ($numItemsDeleted > 0) {
    print STDOUT qq{$pn: Note: the trash directory will be auto-deleted every reboot, without notification.\n};
    print STDOUT (qq{-} x 80) . "\n";
}
