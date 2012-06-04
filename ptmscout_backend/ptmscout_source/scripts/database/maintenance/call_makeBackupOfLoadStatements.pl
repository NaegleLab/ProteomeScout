#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;


my $dbh = returnProductionDBHNOCommit();

my $f = makeBackupOfLoadStatements($dbh);
print "load statements at $f\n";

$dbh->rollback();
