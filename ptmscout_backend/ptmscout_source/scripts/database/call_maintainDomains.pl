#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use DBTools::deleteFunctions;
use errorHandling;
use HMMTools;

# if(scalar(@ARGV) != 1){
#     print "Usage: call_deleteExperiment.pl EXPERIMENT_ID\n";
#     exit;
# }

my $dbh = returnProductionDBHNOCommit();
flushLog();

#my $expId = $ARGV[0];

maintainDomains($dbh);

$dbh->commit();
