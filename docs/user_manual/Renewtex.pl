#!/usr/bin/env perl
while (<STDIN>) {
  $line = $_;
  if ($line =~ /^\\noindent \\textbf/) {
# find first and last curly braces
    $len = length($line);
    $i = 0;
    $ifirst = $len;
    $ilast = -1;
    while ($i < $len) {
      $char = substr($line,$i,1);
      if ($char eq "{") {
        $i++;
        $ifirst = $i;
        last;
      }
      $i++;
    }
    $cnt = 1;
    while ($i < $len) {
      $charm = substr($line,$i-1,1);
      $char = substr($line,$i,1);
      if ($char eq "{" && $charm ne "\\") {
        $cnt++;
      }
      if ($char eq "}" && $charm ne "\\") {
        $cnt--;
      }
      if ($cnt == 0) {
        $i++;
        $ilast = $i;
        last;
      }
      $i++;
    }
    $delta = $ilast-$ifirst-1;
    $tmp = substr($line,$ifirst,$delta);
    $tmp2 = "";
    if (length($line) > $ilast+1) {
      $tmp2 = substr($line,$ilast);
    }
#    $line = "\\textcolor\{red\}\{\\textbf\{\\tt\{$tmp\}\}$tmp2\}\n";
    $line = "\\textcolor\{red\}\{\\texttt\{\\textbf\{$tmp\}\}$tmp2\}\n";
  } elsif ($line =~ /(.*)\\textbf(.*)/) {
#
# Line is not a code block but contains inline code samples
# There may be more than one inline code sample
# Note that @- and @+ are the starting and ending indices of a Perl match
#
    chomp($line);
    $tmp = $line;
    $newline = "\n";
#
# Note that Perl matches last instance of string
#
    while ($tmp =~/(.*)\\textbf\{(.*)/) {
      $tmp = $1;
      $tmp2 = $2;
#
#  Search through tmp2 for closing bracket
#
      $len = length($tmp2);
      $i = 0;
      $ifirst = 0;
      $ilast = -1;
      $cnt = 1;
      while ($i < $len) {
        $charm = "";
        if ($i>0) {
          $charm = substr($tmp2,$i-1,1);
        }
        $char = substr($tmp2,$i,1);
        if ($char eq "{" && $charm ne "\\") {
          $cnt++;
        }
        if ($char eq "}" && $charm ne "\\") {
          $cnt--;
        }
        if ($cnt == 0) {
          $i++;
          $ilast = $i;
          last;
        }
        $i++;
      }
      $tmp3 = "";
      if ($ilast > 0) {
        $tmp3 = substr($tmp2,0,$ilast-1);
        $tmp4 = substr($tmp2,$ilast);
        $newline = "\\texttt\{\\textbf\{"."$tmp3"."\}\}"."$tmp4"."$newline";
      }
    }
    $line = $tmp.$newline;
  }
  if ($line =~ /^\s*\\noindent\s*(.*)/) {
    $line = $1;
  }
  print "$line";
#include color package
  if ($line =~ /usepackage.dvips..graphicx/) {
    print "\\usepackage\{color\}\n";
    print "\\usepackage\{courier\}\n";
    print "\n";
    print "\\setlength\{\\parskip\}\{1ex\}\n";
    print "\\setlength\{\\parindent\}\{0ex\}\n";
  }
}
