#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;

if(scalar(@ARGV) !=2){
    print "USAGE: call_loadMultipleExperimentFile.pl EXPERIMENT_FILE DB\n";
    exit;
}

my($file, $DB) = @ARGV;
my $dbh = returnDBHNOCommit($DB);

loadMultipleExperimentFile($dbh, $file);

$dbh->commit();
