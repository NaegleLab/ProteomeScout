#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use DBTools::deleteFunctions;
use errorHandling;
use HMMTools;

if(scalar(@ARGV) != 2){
    print "Usage: call_deleteExperiment.pl EXPERIMENT_ID DB\n";
    exit;
}

my ($expId, $DB) = @ARGV;

my $dbh = returnDBHNOCommit($DB);
flushLog();


deleteExperiment($dbh, $expId);

$dbh->commit();
