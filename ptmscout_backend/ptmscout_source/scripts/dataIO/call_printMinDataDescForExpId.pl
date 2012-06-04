#! /usr/bin/perl
use strict;
use warnings;
use dataIO;
use DBTools::dbIO;


if(scalar(@ARGV) != 4){

    print "Usage: call_printMinDataDescForExpId.pl EXP_ID OUTPUT_FILE DB_NAME AVERAGE\n";
    exit;
}
my ($expId, $outputFile, $DB, $AVG) = @ARGV; 
 
#my $dbh = returnProductionDBHNOCommit();
my $dbh = returnDBHNOCommit($DB); 

printMinDataDescForExpId($dbh, $expId, $outputFile, $AVG);

$dbh->rollback;
