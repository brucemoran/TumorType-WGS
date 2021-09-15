#! perl

use strict;
use warnings;

##join ../../pcawg_mutation*csv files
open(DIST, "../../pcawg_mutation_distribution.csv");
my %dist;
while(<DIST>){
  chomp $_;
  my @sp=split(/\,/);
  $dist{$sp[0]}=$_;
}
close(DIST)

open(TYPE, "../../pcawg_mutations_types.csv");
open(OUT, "<", "example_input.csv");
while(<TYPE>){
  chomp $_;
  my @sp=split(/\,/);
  print OUT $dist{$sp[0]} . "," . join(",", @sp[1..$#sp]) . "\n";
}
close TYPE;
close OUT;
