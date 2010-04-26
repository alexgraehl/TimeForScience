#!/usr/bin/perl


# Note: if you are including this file, and you want to convert HTML entities like space->%20 (or you want to convert back), then there is a standard perl module that does this:  use URI::Escape;



sub agwUnconditionallyUntaint($) {
    # Insists to perl that a variable is SAFE and not
    # tainted, no matter what. Don't use this without
    # understanding it, or you will have all your files
    # removed by a malicious web user! This is important
    # only when running perl in "taint mode" with -T

    # Info at: http://www.webreference.com/programming/perl/taint/
    if ($_[0] =~ /^(.*)$/) { $_[0] = $1; }
}

sub agwConditionalUntaint($$) {
    # First arg: a variable to un-taint
    # Second: a STRING of something to match. Note: NOT a regular expression! Instead, it is the string that does INTO the regular expression.
    # Here is an example:
    # agwConditionalUntaint($myVar, "A-Za-z\t");
    # Bad example for the second arg: /A-Za-z\t/ <-- NOT like this
    if ($_[0] =~ /^($_[1])$/) {
	$_[0] = $1;
	return 1; # it worked
    } else {
	$_[0] = "TAINTED DATA FOR A VARIABLE was PASSED INTO agwConditionalUntaint. It did not match the regular expression \/\^\($_[1]\)\$\/";
	return 0; # problem
	#warn("TAINTED DATA SENT BY $ENV{'REMOTE_ADDR'}: $drug: $!");
    }

}




sub html2table($$$$)
{
   my ($page, $keep_blanks, $keep_format, $keep_links) = @_;

   my @keep_format_chars  = split(" ", "br b font sup pre");
   my $delim         = "\t";

   $page =~ s/<![^>]*>//gi;

   $page =~ s/<\/html[^>]*>//gi;

   $page =~ s/<\/body[^>]*>//gi;

   my @lines = split(/\<table[^>]*>/i, $page);

   my @tables;

   for(my $t = 1; $t < @lines; $t++)
   {
      my @table;

      $lines[$t] =~ s/<\/table[^>]*>//gi;

      my @rows = split(/<tr[^>]*>/i, $lines[$t]);

      foreach my $row (@rows)
      {
         if($keep_blanks or $row =~ /<td/i)
         {
            $row =~ s/<td[^>]*>/$delim/gi;
            $row =~ s/\n//g;
            $row =~ s/<\/td[^>]*>//gi;
            $row =~ s/<\/tr[^>]*>//gi;

            if(not($keep_format))
            {
               foreach my $keep_format_char (@keep_format_chars)
               {
                  $row =~ s/<$keep_format_char[^>]*>/ /gi;
                  $row =~ s/<\/$keep_format_char[^>]*>//gi;
               }
            }

            if(not($keep_links))
            {
               # Get rid of links
               $row =~ s/<a\s+href\s*=[^>]*>([^<]*)<\/a[^>]*>/$1/gi;

               # Get rid of images
               $row =~ s/<\s*img[^>]*>//gi;

               # Get rid of labels
               $row =~ s/<a\s+Name\s*=[^>]*>//gi;
            }
            if($keep_blanks or $row =~ /\S/)
            {
               my @row = split($delim, $row);
               push(@table, \@row);
            }
         }
      }
      if($keep_blanks or scalar(@table) > 0)
      {
         push(@tables, \@table);
      }
   }
   return \@tables;
}


#-------------------------------------------------------------------------------
# $string getHtmlTableRow(\@list ; $color)
#-------------------------------------------------------------------------------
sub getHtmlTableRow(\@;$) {
    # color is an OPTIONAL second argument for the cell color for each cell on this row

   my ($list,$color) = @_;
   $list = defined($list) ? $list : [split];

   my $colorStr = '';
   if (defined($color) and $color) {
       $colorStr = qq{ BGCOLOR=\"$color\"};
   }

   my $string = "<TR${colorStr}>\n";

   foreach my $entry (@{$list}) {
      $string .= "     <TD>${entry}</TD>\n";
   }
   $string .= "</TR>";

   return $string;
}


#-------------------------------------------------------------------------------
# $string getHtmlTableRow(\@list ; $rowAttributes $colAttributes)
#-------------------------------------------------------------------------------
sub getFormattedHtmlTableRow(\@;$$) {
    # $rowAttr: Row attributes: something that is done to the TR tag (text that goes in the <TR ... >
    # $cellAttr: Cell attributes: something that is done to each cell (text that goes in the <TD ... >

   my ($list,$rowAttr,$cellAttr) = @_;
   $list = defined($list) ? $list : [split];

   if (!defined($rowAttr) || (scalar($rowAttr)>0)) { $rowAttr = ''; } else { $rowAttr = ' ' . $rowAttr; } # put a space before the attribute, if there IS one defined
   if (!defined($cellAttr) || (scalar($cellAttr)>0)) { $cellAttr = ''; } else { $cellAttr = ' ' . $cellAttr; } # put a space before the attribute, if there IS one defined

   my $string = "<TR${rowAttr}>\n"; # <-- the space gets added above so that there is a space after the TR

   foreach my $entry (@{$list}) {
      $string .= "     <TD${cellAttr}>${entry}</TD>\n";
   }
   $string .= "</TR>";

   return $string;
}






# =========== end of content (1 afterward to make Perl not freak out when you include it)
1  # <-- keep the 1 here!
