#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbEntry;
use DBTools::dbIO;

if (scalar(@ARGV) < 1){
    print "Usage: call_handleAllProteins.pl ACC1 ACC2.. ACCN\n";
    print "INTO TEST DATABASE\n";
    exit;

}

my $dbh = returnTestDBH();
my @acc = @ARGV;
my($proteinId, $iProteinFlag, $iAccFlag);

my ($hashRef) = handleAllProteins($dbh, \@acc);
#for(my $i=0; $i<=scalar(@acc); $i++){
 foreach my $acc (keys %$hashRef){
     my $arr = $hashRef->{$acc};
    print "Accession: $acc\tProteinId: $arr->[0]\tiProteinFlag: $arr->[1]\tiAccFlag: $arr->[2]\n";

}

