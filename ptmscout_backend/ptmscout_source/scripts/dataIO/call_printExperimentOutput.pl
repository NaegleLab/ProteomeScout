#! /usr/bin/perl
use strict;
use warnings;
use dataIO;
use DBTools::queryTools;
use DBTools::dbIO;

if(scalar(@ARGV) != 2){
    print "Usage: call_printExperimentOutput.pl OUTPUT_FILE EXPERIMENT_ID\n";
    exit;

}

my $dbh = returnProductionDBHNOCommit();
my ($outputFile, $expId) = @ARGV;

#my $outputFile = 'test_printExperiment_out.txt';
#my $expId = 22;
printExperimentOutput($dbh, $outputFile, $expId);

$dbh->rollback();
