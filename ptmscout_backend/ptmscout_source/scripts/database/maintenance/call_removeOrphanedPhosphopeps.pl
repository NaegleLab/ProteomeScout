#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;

my $dbh = returnProductionDBHNOCommit();

my $phosphopeps = removeOrphanedPhosphopeps($dbh);
print "Found ".scalar(@$phosphopeps)." orphaned peptides\n";

$dbh->commit();
