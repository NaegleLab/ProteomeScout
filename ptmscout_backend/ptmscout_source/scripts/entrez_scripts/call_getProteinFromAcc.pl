#! /usr/bin/perl
use strict;
use warnings;
use entrezTools;


if(scalar(@ARGV) < 1){
    print "Usage: call_getProteinFromAcc.pl GI\n";
    exit;
}
my $gi = $ARGV[0];
my ($errorCode, $sequence, $species, $gene, $name, $locus);


($errorCode, $sequence, $species, $gene, $name, $locus) = getProteinFromAcc($gi);
print "$errorCode\n $sequence\n $species\n $gene\n $name\n $locus\n";
