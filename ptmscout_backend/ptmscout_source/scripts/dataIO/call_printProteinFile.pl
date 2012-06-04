#! /usr/bin/perl
use strict;
use warnings;
use dataIO;
use DBTools::dbIO;


if(scalar(@ARGV) != 2){
    print "USAGE: call_printProteinFile.pl EXP_ID OUTPUT_FILE\n";
    exit;
}

my ($expId, $outputFile) = @ARGV;

my $dbh = returnProductionDBHNOCommit();
printProteinFile($dbh, $expId, $outputFile);

$dbh->rollback();
