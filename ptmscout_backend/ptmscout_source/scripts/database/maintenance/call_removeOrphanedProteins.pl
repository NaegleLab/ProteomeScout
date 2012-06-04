#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::queryTools;
use DBTools::maintenance;

my $dbh = returnProductionDBHNOCommit();


my $removed = removeOrphanedProteins($dbh);
print "Removed ".scalar(@$removed)." proteins from database\n";

$dbh->commit();
