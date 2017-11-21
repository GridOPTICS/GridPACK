#!/usr/bin/env perl
#
# Start parsing file. Assume that top of file is boiler plate and ignore
# it until you hit the first preprocessor statement.
#
$ignore = 1;
$isComment = 0;
$isPreprocess = 0;
$isDox = 0;
$isCode = 0;
$firstCodeLine = 0;
$firstDoxLine = 0;
$firstLatexLine = 0;
print "\\documentclass\[12pt\]\{article\}\n";
print "\\usepackage\{color\}\n";
print "\\usepackage\{alltt\}\n";
print "\\usepackage\{lmodern\}\n";
print "\\begin\{document\}\n";
print "\\setlength\{\\oddsidemargin\}\{-1.0cm\}\n";
print "\\setlength\{\\textwidth\}\{17cm\}\n";

while (<STDIN>) {
  $line = $_;
  if (/^\#/) {
    $ignore = 0;
  }
  if ($ignore == 0) {
    if (/^\s*\#/) {
#
# C++/C preprocessor lines (starts with #)
#
      $isPreprocess = 1;
      $isComment = 0;
      $isDox = 0;
      if ($isCode) {
        print "\\end\{alltt\}\n";
      }
      $isCode = 0;
      $firstCodeLine = 1;
      $firstDoxLine = 1;
      $firstLatexLine = 1;
    } elsif (/^\s*\/\*\*/) {
#
# Doxygen comment lines lines (starts with /**)
#
      $isPreprocess = 0;
      $isComment = 0;
      $isDox = 1;
      if ($isCode) {
        print "\\end\{alltt\}\n";
      }
      $isCode = 0;
      $firstCodeLine = 1;
      $firstLatexLine = 1;
    } elsif (/\<latex\>/) {
#
# Inline LaTex documentation (starts with <latex>) (behaves
# similarly to Doxygen documentation)
#
      $isPreprocess = 0;
      $isComment = 1;
      $isDox = 0;
      if ($isCode) {
        print "\\end\{alltt\}\n";
      }
      $isCode = 0;
      $firstDoxLine = 1;
      $firstCodeLine = 1;
    } elsif (!$isPreprocess && !$isComment && !$isDox) {
      $isPreprocess = 0;
      $isComment = 0;
      $isDox = 0;
      $isCode = 1;
      $firstDoxLine = 1;
      $firstLatexLine = 1;
    }
    chomp($line);
    if ($isPreprocess) {
      if ($line =~ /(^\#)(\s*\S+)(.*)/) {
        $line =  "\\textcolor\{magenta\}\{$1$2\}\\textcolor\{red\}\{$3\}";
      }
      print "\\begin\{alltt\}\n";
      print "$line\n";
      print "\\end\{alltt\}\n";
      $isPreprocess = 0;
    }
    if ($isComment) {
      if ($firstLatexLine) {
        printf "\\texttt\{\\textcolor\{cyan\}\{";
        $firstLatexLine = 0;
      }
      # Remove opening delimiter
      if ($line =~ /\<latex\>\s*(\S.*)/) {
        $line = $1;
      }
      # Look for closing delimiter
      if ($line =~ /(.*)\<\/latex\>/) {
        $isComment = 0;
        $line = "$1\}\}";
      }
      # Remove C++ comment delimiters
      if ($line =~ /\s*\/\/\s*(\S.*)/) {
        $line = $1;
      }
      print "$line\n";
    }
    if ($isCode) {
      if ($firstCodeLine) {
         print "\\begin\{alltt\}\n";
         $firstCodeLine = 0;
      }
      # Search for comments
      if ($line =~ /(.*)(\/\/)(.*)/) {
        $one = $1;
        $two = $2;
        $three = $3;
        # Fix up some symbols that cause problems
        $one =~ s/\}/\\\}/g;
        $one =~ s/\\/\\textbackslash\\hspace\{0pt\}/g;
        $one =~ s/\{/\\\{/g;
        $two =~ s/\\/\\textbackslash\\hspace\{0pt\}/g;
        $two =~ s/\{/\\\{/g;
        $two =~ s/\}/\\\}/g;
        $three =~ s/\\/\\textbackslash\\hspace\{0pt\}/g;
        $three =~ s/\{/\\\{/g;
        $three =~ s/\}/\\\}/g;
        $line = "$one\\textcolor\{blue\}\{$two$three\}";
      } else {
        $line =~ s/\\/\\textbackslash\\hspace\{0pt\}/g;
        # Fix up some symbols that cause problems
        $line =~ s/\{/\\\{/g;
        $line =~ s/\}/\\\}/g;
      }
     $line =~ s/\\\{0pt\\\}/\{0pt\}/g;
      print "$line\n";
    }
    if ($isDox) {
      if ($firstDoxLine) {
        printf "\\texttt\{\\textcolor\{cyan\}\{";
        $firstDoxLine = 0;
      }
      # Look for closing delimiter
      if ($line =~ /(.*)\*\//) {
        $isDox = 0;
        $line = $1;
	if ($line =~/(\*+)([^\*]*)/) {
	  $line = $2;
	}
	if ($line =~/([^\*]*)(\*+)/) {
	  $line = $1;
	}
        if (!($line =~ /\s*/)) {
          $line = "$chk\}\}\n";
        } else {
          $line = "\}\}";
        }
        if ($line =~ /\@param\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@param\} $1\} $2";
        }
        print "\\noindent\n";
        print "$line\n";
      } else {
#
#  Get rid of any remaining asterisks
#
	if ($line =~/(\*+)([^\*]*)/) {
	  $line = $2;
	}
	if ($line =~/([^\*]*)(\*+)/) {
	  $line = $1;
	}
        if ($line =~ /\@param\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@param\} $1\} $2";
        }
        print "\\noindent\n";
        print "$line\\\\\n";
      }
    }
    if ($isLatex) {
      if ($firstDoxLine) {
        printf "\\texttt\{\\textcolor\{cyan\}\{";
        $firstDoxLine = 0;
      }
      # Look for closing delimiter
      if ($line =~ /(.*)\*\//) {
        $isDox = 0;
        $line = $1;
	if ($line =~/(\*+)([^\*]*)/) {
	  $line = $2;
	}
	if ($line =~/([^\*]*)(\*+)/) {
	  $line = $1;
	}
        if (!($line =~ /\s*/)) {
          $line = "$chk\}\}\n";
        } else {
          $line = "\}\}";
        }
        if ($line =~ /\@param\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@param\} $1\} $2";
        }
        print "\\noindent\n";
        print "$line\n";
      } else {
#
#  Get rid of any remaining asterisks
#
	if ($line =~/(\*+)([^\*]*)/) {
	  $line = $2;
	}
	if ($line =~/([^\*]*)(\*+)/) {
	  $line = $1;
	}
        if ($line =~ /\@param\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@param\} $1\} $2";
        }
        print "\\noindent\n";
        print "$line\\\\\n";
      }
    }
  }
}
if ($isCode) {
  print "\\end\{alltt\}\n";
}
print "\\end\{document\}\n";
