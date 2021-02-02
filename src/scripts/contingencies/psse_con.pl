#!/usr/bin/env perl
#
#  This script assumes that data in the .raw file is in a v33 format.
#
#  Read in command line arguments. This script requires at least 2 arguments and
#  optionally several more. The two required arguments are the name of the PSS/E
#  .RAW file that is used to generate the contingencies and cost parameter files
#  and the second argument is the root name that will be used for the output
#  files.
#
#  The optional directives are
#
#  -GC  (Single generator contingencies)
#  -LC  (Single line contingencies)
#  -GGC (Two generator contingencies)
#  -GLC (A generator and a line contingecies)
#  -LLC (Two line contingencies)
#
#  Each of these directives can be used to control the number of contingencies
#  of a given type that are added to the contingency file.
#
#  The allowable arguments for these directives are 
#  
#  all
#  none
#  some number
#
#  If "all is specified, then all possible contingencies of that type will be
#  generated, if "none" is specified then no contingencies of that type will be
#  generated and if a number is specified, then that number of contingencies of
#  that type is generated. By default, the -GC and -LC options are set to "all"
#  and the -GGC, -GLC and -LLC options are set to none
#
$argc = @ARGV;
$rawfile = "";
if ($argc > 0) {
  $rawfile = $ARGV[0];
} else {
  print("No RAE file specified. Exiting...\n");
  exit(0);
}
print "PSS/E FILE NAME: $rawfile\n";
$rootname = "";
if ($argc > 1) {
  $rootname = $ARGV[1];
}
print "EXPORT ROOT: $rootname\n";
#
#  Parse arguments to find out how many contingencies to create. If no arguments
#  create all N-1 generator and line contingencies and no N-2 contingencies.
#
$allGC = 1;
$allLC = 1;
$allGGC = 0;
$allGLC = 0;
$allLLC = 0;
$maxGC = 0;
$maxLC = 0;
$maxGGC = 0;
$maxGLC = 0;
$maxLLC = 0;
for ($i=0; $i<$argc; $i++) {
  if ($ARGV[$i] eq "-GC" && defined($ARGV[$i+1])) {
    $arg = lc($ARGV[$i+1]);
    if ($arg eq "all") {
      $allGC = 1;
    } elsif ($arg eq "none") {
      $allGC = 0;
    } else {
      $allGC = 0;
      $maxGC = $arg;
    }
  } elsif ($ARGV[$i] eq "-LC" && defined($ARGV[$i+1])) {
    $arg = lc($ARGV[$i+1]);
    if ($arg eq "all") {
      $allLC = 1;
    } elsif ($arg eq "none") {
      $allLC = 0;
    } else {
      $allLC = 0;
      $maxLC = $arg;
    }
  } elsif ($ARGV[$i] eq "-GGC" && defined($ARGV[$i+1])) {
    $arg = lc($ARGV[$i+1]);
    if ($arg eq "all") {
      $allGGC = 1;
    } elsif ($arg eq "none") {
      $allGGC = 0;
    } else {
      $allGGC = 0;
      $maxGGC = $arg;
    }
  } elsif ($ARGV[$i] eq "-GLC" && defined($ARGV[$i+1])) {
    $arg = lc($ARGV[$i+1]);
    if ($arg eq "all") {
      $allGLC = 1;
    } elsif ($arg eq "none") {
      $allGLC = 0;
    } else {
      $allGLC = 0;
      $maxGLC = $arg;
    }
  } elsif ($ARGV[$i] eq "-LLC" && defined($ARGV[$i+1])) {
    $arg = lc($ARGV[$i+1]);
    if ($arg eq "all") {
      $allLLC = 1;
    } elsif ($arg eq "none") {
      $allLLC = 0;
    } else {
      $allLLC = 0;
      $maxLLC = $arg;
    }
  }
}
# 
$count = 1;
#
#  Open file for contingencies
#
$filename = "$rootname\.con";
#
#  Parse PSSE format file for information on generators and lines
#
open(RAW,"$rawfile");
open(CONTINGENCY,">$filename");
print CONTINGENCY "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
print CONTINGENCY "<ContingencyList>\n";
print CONTINGENCY "  <Contingencies>\n";
$is_gen = 0;
$is_branch = 0;
$is_bus = 0;
$is_xform = 0;
$is_gen_parsed = 0;
$is_branch_parsed = 0;
$is_bus_parsed = 0;
@nminus2 = ();
$count2 = 0;
$gcnt = 0;
$lcnt = 0;
$count_gc = 0;
$count_lc = 0;
$buscnt = 0;
@fromBus = ();
@toBus = ();
@lineTags = ();
%linecnt = {};
%busmap = {};
@busIDs = ();
# read first three lines to get to bus block
$line = <RAW>;
$line = <RAW>;
$line = <RAW>;
$is_bus = 1;
while (<RAW>) {
  $line = $_;
  if ($line =~ /END OF BUS DATA/ && $is_bus == 1) {
    $is_bus = 0;
    $is_bus_parsed = 1;
    print "Closing buses: $line";
  }
  if ($line =~ /BEGIN GENERATOR DATA/) {
    print "Begin generator data\n";
    $is_gen = 1;
    $line = <RAW>;
  }
  if ($line =~ /END OF GENERATOR DATA/ && $is_gen == 1) {
    $is_gen = 0;
    $is_gen_parsed = 1;
    print "Closing generators: $line";
  }
  if ($line =~ /BEGIN BRANCH DATA/) {
    print "Begin branch data\n";
    $is_branch = 1;
    $line = <RAW>;
  }
  if ($line =~ /END OF BRANCH DATA/ && $is_branch == 1) {
    $is_branch = 0;
    print "Closing lines: $line";
  }
  if ($line =~ /BEGIN TRANSFORMER DATA/) {
    $is_xform = 1;
    $line = <RAW>;
  }
  if ($line =~ /END OF TRANSFORMER DATA/ && $is_xform == 1) {
    $is_xform = 0;
    $is_branch_parsed = 1;
    print "Closing lines: $line";
  }
  if ($is_gen == 1) {
    @tokens = split(',', $line);
    $busid = $tokens[0];
    $busid =~ s/^\s*//;
    $busid =~ s/\s*$//;
    $tag = $tokens[1];
    $tag =~ s/^\s*\'//;
    $tag=~ s/\s*\'\s*$//;
    $active = $tokens[14];
    $active =~ s/^\s*//;
    $active =~ s/\s*$//;
    if ($active > 0) {
      if ($allGC == 1 || ($count_gc < $maxGC)) {
        print CONTINGENCY "    <Contingency>\n";
        print CONTINGENCY "      <contingencyType> Generator <\/contingencyType>\n";
        print CONTINGENCY "      <contingencyName> CTG$count <\/contingencyName>\n";
        print CONTINGENCY "      <contingencyBuses> $busid <\/contingencyBuses>\n";
        print CONTINGENCY "      <contingencyGenerators> $tag <\/contingencyGenerators>\n";
        print CONTINGENCY "    <\/Contingency>\n";
        $contingency = "GENERATOR\_$busid\_$tag";
        $nminus2[$count2] = $contingency;
        $count++;
        $count_gc++;
        $count2++;
      }
      $gcnt++;
    }
  }
#
# Just gather line information for this pass. Construct contingencies on second
# pass
#
  if ($is_branch == 1) {
    @tokens = split(',',$line);
    $bus1 = $tokens[0];
    $bus1 =~ s/^\s*//;
    $bus1 =~ s/\s*$//;
    $bus2 = $tokens[1];
    $bus2 =~ s/^\s*//;
    $bus2 =~ s/\s*$//;
    $tag = $tokens[2];
    $tag =~ s/^\s*\'//;
    $tag =~ s/\s*\'\s*$//;
    $active = $tokens[13];
    $active =~ s/^\s*//;
    $active =~ s/\s*$//;
    if ($active > 0) {
      if (defined($linecnt{$bus1})) {
        $linecnt{$bus1}++;
      } else {
        $linecnt{$bus1} = 1;
      }
      if (defined($linecnt{$bus2})) {
        $linecnt{$bus2}++;
      } else {
        $linecnt{$bus2} = 1;
      }
      $fromBus[$lcnt] = $bus1;
      $toBus[$lcnt] = $bus2;
      $lineTags[$lcnt] = $tag;
      $lcnt++;
    }
  }
  if ($is_xform == 1) {
    @tokens = split(',',$line);
    $bus1 = $tokens[0];
    $bus1 =~ s/^\s*//;
    $bus1 =~ s/\s*$//;
    $bus2 = $tokens[1];
    $bus2 =~ s/^\s*//;
    $bus2 =~ s/\s*$//;
    $flag = $tokens[2];
    $flag =~ s/^\s*//;
    $tflg =~ s/\s*$//;
    $tag = $tokens[3];
    $tag =~ s/^\s*//;
    $tag =~ s/\s*$//;
    $active = $tokens[11];
    $active =~ s/^\s*//;
    $active =~ s/\s*$//;
    if ($active > 0) {
      if (defined($linecnt{$bus1})) {
        $linecnt{$bus1}++;
      } else {
        $linecnt{$bus1} = 1;
      }
      if (defined($linecnt{$bus2})) {
        $linecnt{$bus2}++;
      } else {
        $linecnt{$bus2} = 1;
      }
      $fromBus[$lcnt] = $bus1;
      $toBus[$lcnt] = $bus2;
      $lineTags[$lcnt] = $tag;
      $lcnt++;
    }
    if ($flag == 0) {
      $line = <RAW>;
      $line = <RAW>;
      $line = <RAW>;
    } else {
      $line = <RAW>;
      $line = <RAW>;
      $line = <RAW>;
      $line = <RAW>;
    }
  }
  if ($is_bus == 1) {
    @tokens = split(',',$line);
    $bus = $tokens[0];
    $bus =~ s/^\s*//;
    $bus =~ s/\s*$//;
    $busmap{$bus} = $buscnt;
    $busIDs[$buscnt] = $bus;
    $buscnt++;
  }
}
print "Buses: $buscnt Generators: $gcnt Active Lines: $lcnt\n";
#
#  create linked list of neighbors
#
@top = ();
@next = ();
@nghbrs = ();
@nghbrtag = ();
for ($ibus = 0; $ibus<$buscnt; $ibus++) {
  $top[$ibus] = -1;
}
$ncnt = 0;
for ($iline = 0; $iline<$lcnt; $iline++) {
  $bus1 = $fromBus[$iline];
  $bus2 = $toBus[$iline];
  $fbus = $busmap{$bus1};
  $tbus = $busmap{$bus2};
  $tmp = $top[$fbus];
  $top[$fbus] = $ncnt;
  $next[$ncnt] = $tmp;
  $nghbrs[$ncnt] = $tbus;
  $nghbrtag[$ncnt] = $lineTags[$iline];
  $ncnt++;
  $tmp = $top[$tbus];
  $top[$tbus] = $ncnt;
  $next[$ncnt] = $tmp;
  $nghbrs[$ncnt] = $fbus;
  $nghbrtag[$ncnt] = $lineTags[$iline];
  $ncnt++;
}
print "Total neighbor count: $ncnt\n";
#
#  Set up line contingencies
#
for ($iline = 0; $iline < $lcnt; $iline++) {
  $bus1 = $fromBus[$iline];
  $bus2 = $toBus[$iline];
  $tag = $lineTags[$iline];
  $status = 1;
  if (!defined($linecnt{$bus1}) || $linecnt{$bus1} < 2) {
    $status = 0;
  }
  if (!defined($linecnt{$bus2}) || $linecnt{$bus2} < 2) {
    $status = 0;
  }
  if ($status == 1) {
    if ($allGC == 1 || ($count_lc < $maxLC)) {
#
#  Check to see if contingency splits network into two pieces. Loop over
#  neighbors recursively until you run out of new neighbors. If total number of
#  buses in cluster is less than the total number of buses in the original
#  network, go to next line outage
#
      $bcnt = 0;
      @newbus = ();
      @found = ();
      for ($ibus = 0; $ibus < $buscnt; $ibus++) {
        $found[$ibus] = 0;
      }
#
#  Find first bus
#
      $ifirst = 0;
      $fbus = $busmap{$bus1};
      $tbus = $busmap{$bus2};
      while ($ifirst+1 == $fbus || $ifirst+1 == $tbus) {
        $ifirst++;
      }
      $found[$ifirst] = 1;
      $bcnt = 1;
      $nbus = 0;
      $ibus = $top[$ifirst];
      while ($ibus >= 0) {
        $j = $nghbrs[$ibus];
        $newbus[$nbus] = $j;
        $found[$j] = 1;
        $bcnt++;
        $nbus++;
        $ibus = $next[$ibus];
      }
      while ($nbus > 0) {
        @oldbus = ();
        for ($i=0; $i<$nbus; $i++) {$oldbus[$i] = $newbus[$i];}
        $newcnt = 0;
        @newbus = ();
        for ($i=0; $i<$nbus; $i++) {
          $j = $top[$oldbus[$i]];
          while ($j >= 0) {
            if (!(($oldbus[$i]+1 == $fbus && $nghbrs[$j]+1 == $tbus && $nghbrtag[$j] == $tag)
              || ($oldbus[$i]+1 == $tbus && $nghbrs[$j]+1 == $fbus && $nghbrtag[$j] == $tag))) {
              $ob = $oldbus[$i]+1;
              $nb = $nghbrs[$j]+1;
              if ($found[$nghbrs[$j]] == 0) {
                $bcnt++;
                $found[$nghbrs[$j]] = 1;
                $newbus[$newcnt] = $nghbrs[$j];
                $newcnt++;
              }
            }
            $j = $next[$j];
          }
        }
        $nbus = $newcnt;
      }
      if ($bcnt == $buscnt) {
        print CONTINGENCY "    <Contingency>\n";
        print CONTINGENCY "      <contingencyType> Line <\/contingencyType>\n";
        print CONTINGENCY "      <contingencyName> CTG$count <\/contingencyName>\n";
        print CONTINGENCY "      <contingencyLineBuses> $bus1 $bus2 <\/contingencyLineBuses>\n";
        print CONTINGENCY "      <contingencyLineNames> $tag <\/contingencyLineNames>\n";
        print CONTINGENCY "    <\/Contingency>\n";
        $contingency = "LINE\_$bus1\_$bus2\_$tag";
        $nminus2[$count2] = $contingency;
        $count_lc++;
        $count2++;
        $count++;
      } else {
        print "Rejecting contingency $bus1-$bus2 cnt: $bcnt nbus: $buscnt\n";
      }
    }
  }
}

#
# print out N-2 contingencies
#
$ncont = @nminus2;
$ngg = 0;
$ngl = 0;
$nll = 0;
for ($i=0; $i<$ncont; $i++) {
  for ($j=$i+1; $j<$ncont; $j++) {
    $type = "";
    @buses = ();
    @gtags = ();
    @ltags = ();
    $nbus = 0;
    $ngtag = 0;
    $nltag = 0;
    if ($nminus2[$i] =~/GENERATOR/) {
      $type = "G";
      if ($nminus2[$i] =~/GENERATOR\_(.*)\_(.*)/) {
         $buses[$nbus] = $1;
         $nbus++;
         $gtags[$ngtag] = $2;
         $ngtag++;
      }
    } elsif ($nminus2[$i] =~/LINE/) {
      $type = "L";
      if ($nminus2[$i] =~/LINE\_(.*)\_(.*)\_(.*)/) {
         $buses[$nbus] = $1;
         $nbus++;
         $buses[$nbus] = $2;
         $nbus++;
         $ltags[$nltag] = $3;
         $nltag++;
      }
    }
    if ($nminus2[$j] =~/GENERATOR/) {
      $type .= "G";
      if ($nminus2[$j] =~/GENERATOR\_(.*)\_(.*)/) {
         $buses[$nbus] = $1;
         $nbus++;
         $gtags[$ngtag] = $2;
         $ngtag++;
      }
    } elsif ($nminus2[$j] =~/LINE/) {
      $type .= "L";
      $contingency = "LINE\_$bus1\_$bus2\_$tag";
      if ($nminus2[$j] =~/LINE\_(.*)\_(.*)\_(.*)/) {
         $buses[$nbus] = $1;
         $nbus++;
         $buses[$nbus] = $2;
         $nbus++;
         $ltags[$nltag] = $3;
         $nltag++;
      }
    }
    if ($type eq "GG") {
      if ($allGGC == 1 || $ngg < $maxGGC) {
        print CONTINGENCY "    <Contingency>\n";
        print CONTINGENCY "      <contingencyType> Generator-Generator <\/contingencyType>\n";
        print CONTINGENCY "      <contingencyName> CTG$count <\/contingencyName>\n";
        print CONTINGENCY "      <contingencyBuses> $buses[0] $buses[1] <\/contingencyBuses>\n";
        print CONTINGENCY "      <contingencyGenerators> $gtags[0] $gtags[1] <\/contingencyGenerators>\n";
        print CONTINGENCY "    <\/Contingency>\n";
        $count++;
        $ngg++;
      }
    } elsif ($type eq "GL") {
      if ($allGLC == 1 || $ngl < $maxGLC) {
        print CONTINGENCY "    <Contingency>\n";
        print CONTINGENCY "      <contingencyType> Generator-Line <\/contingencyType>\n";
        print CONTINGENCY "      <contingencyName> CTG$count <\/contingencyName>\n";
        print CONTINGENCY "      <contingencyBuses> $buses[0] <\/contingencyBuses>\n";
        print CONTINGENCY "      <contingencyGenerators> $gtags[0] <\/contingencyGenerators>\n";
        print CONTINGENCY "      <contingencyLineBuses> $buses[1] $buses[2] <\/contingencyLineBuses>\n";
        print CONTINGENCY "      <contingencyLineNames> $ltags[0] <\/contingencyLineNames>\n";
        print CONTINGENCY "    <\/Contingency>\n";
        $count++;
        $ngl++;
      }
    } elsif ($type eq "LL") {
      if ($allLLC == 1 || $nll < $maxLLC) {
        print CONTINGENCY "    <Contingency>\n";
        print CONTINGENCY "      <contingencyType> Line-Line <\/contingencyType>\n";
        print CONTINGENCY "      <contingencyName> CTG$count <\/contingencyName>\n";
        print CONTINGENCY "      <contingencyLineBuses> $buses[0] $buses[1] $buses[2] $buses[3] <\/contingencyLineBuses>\n";
        print CONTINGENCY "      <contingencyLineNames> $ltags[0]  $ltags[1] <\/contingencyLineNames>\n";
        print CONTINGENCY "    <\/Contingency>\n";
        $count++;
        $nll++;
      }
    }
  }
}
print CONTINGENCY "  <\/Contingencies>\n";
print CONTINGENCY "<\/ContingencyList>\n";
close(CONTINGENCY);
