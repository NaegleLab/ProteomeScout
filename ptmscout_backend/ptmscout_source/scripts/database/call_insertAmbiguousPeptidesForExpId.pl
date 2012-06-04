#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::queryTools;
use DBTools::insertFunctions;
use DBTools::dbEntry;

if(scalar(@ARGV) != 1){
    print "Usage: call_insertAmbiguousPeptidesForExpId.pl EXP_ID\n";
    exit;
}

my $dbh = returnProductionDBHNOCommit();
my $rdbh = returnRefSeqDBHNOCommit();

my $expId = $ARGV[0];
my $count;

#$expId = 8; # MRM dataset
$count = insertAmbiguousPeptidesForExpId($dbh, $rdbh, $expId);
print "There were $count ambiguous peptides\n";



$dbh->commit();
$rdbh->rollback();
