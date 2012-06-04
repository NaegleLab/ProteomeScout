#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::deleteFunctions;

if(scalar(@ARGV) != 1){

    print "Usage: call_deleteProteinRecordsByProteinId.pl PROTEIN_ID\n";
    exit;
}

my ($proteinId) = @ARGV;
my $dbh = returnProductionDBHNOCommit();
my $error = deleteProteinRecordsByProteinId($dbh, $proteinId);
print "ERROR: $error\n";

$dbh->commit();
