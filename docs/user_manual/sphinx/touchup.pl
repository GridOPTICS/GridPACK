#!/usr/bin/env perl
while(<STDIN>) {
  $line = $_;
  $line =~ s/\*\*\`\`/\`\`/g;
  $line =~ s/\`\`\*\*/\`\`/g;
  print "$line";
}
