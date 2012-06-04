#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;
use DBTools::dbEntry;

if(scalar(@ARGV) != 0){

    print "Usage: call_maintainDomains\n";
    exit;
}

my $dbh = returnProductionDBHNOCommit();

maintainDomains($dbh);


$dbh->commit();
