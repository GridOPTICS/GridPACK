#!/usr/bin/env perl
#
# This script is designed to turn C++ files (.cpp and .hpp) into LaTeX formatted
# files that can be used for creating PDFs of the code. The script ignores all
# text before the first preprocessor symbol (a line starting with the # sign) so
# any file that does not have a #include line or some other preprocessor symbol
# at the top will not work.
#
# The script can also be used to include LaTeX equations in the comments. If the
# script contains Doxygen-style comments (beginning with the symbol /**) then
# LaTeX equations can be included in the comment using the $ delimiters. Inline
# comments can also be used to include LaTeX equations. Comments that start with
# <latex> and end with </latex> will be interpreted as LaTeX text. Thus, you can
# include something like
#
# // <latex> These lines evaluate the function $e^{-x^2}sin(x)$ for $x$ in the
# // interval $\[0,X_{max}\]$ </latex>
#
# and get the expected formulas in the final PDF. This will not affect Doxygen
# processing, although the code fragments will not come out correctly. If you
# include <latex> delimited text in in-code comments, these will get processed
# differently from code comments delimited by the usual // symbols. It will
# probably look better if you include math style comments at the top of the
# routine instead of inside it.
#
###############################################################################
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
# Inline LaTeX documentation (starts with <latex>) (behaves
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
      #
      # Convert any Doxygen \f$ signs to just $s
      #
      $line =~ s/\\f\$/\$/g;
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
        $one =~ s/\\/\\textbackslash\\hspace\{0pt\}/g;
        $one =~ s/\}/\\\}/g;
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
        print "\\noindent\n";
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
        if ($line =~ /\@return\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@return\} $1\} $2";
        }
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
        if ($line =~ /\@return\s+(\S+)\s+(.*)/) {
          $line = "\\textbf\{\\textcolor\{red\}\{\@return\} $1\} $2";
        }
        print "$line\\\\\n";
      }
    }
    if ($isLatex) {
      if ($firstDoxLine) {
        printf "\\texttt\{\\textcolor\{cyan\}\{\n";
        print "\\noindent\n";
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
      }
      print "$line\n";
    }
  }
}
if ($isCode) {
  print "\\end\{alltt\}\n";
}
print "\\end\{document\}\n";
