#! /usr/bin/perl
use strict;
use warnings;
use DBI;
use DBTools::queryTools;
use DBTools::dbIO;
use DBTools::insertFunctions;


if(scalar(@ARGV) != 1){
    print "Usage: call_returnGeneListForExperiment.pl EXPERIMENT_ID\n";
    exit;
}

my $dbh = returnProductionDBHNOCommit();
my $expId = $ARGV[0];

my $accHash = returnGeneListForExperiment($dbh, $expId);
print "Protein ID\tAccession\n";
foreach my $protein (keys %$accHash){
    print "$protein\t$accHash->{$protein}\n";
}


$dbh->rollback();
