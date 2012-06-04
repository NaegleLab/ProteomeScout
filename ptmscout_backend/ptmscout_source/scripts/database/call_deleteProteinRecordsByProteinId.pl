#! /usr/bin/perl
use strict;
use warnings;
use DBTools::deleteFunctions;
use DBTools::dbIO;


if(scalar(@ARGV) != 1){
    print "USAGE: call_deleteProteinRecordsByProteinId.pl PROTEIN_ID\n";
    exit;
}

my ($proteinId) = @ARGV;
my $dbh = returnProductionDBHNOCommit();

deleteProteinRecordsByProteinId($dbh, $proteinId);

$dbh->commit();
