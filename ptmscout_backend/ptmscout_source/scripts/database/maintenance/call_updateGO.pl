#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;
use GOTools;

print "\n\nYOU MUST MANUALLY UPDATE the OBO file\n\n\n";

if(scalar(@ARGV) != 1){
    print "USAGE: call_updateGO.pl REMOVE_IEA\n";
    exit;
}

my ($REMOVE_IEA) = @ARGV;

my $dbh = returnProductionDBHNOCommit();

updateGO($dbh, $REMOVE_IEA);
$dbh->commit();
