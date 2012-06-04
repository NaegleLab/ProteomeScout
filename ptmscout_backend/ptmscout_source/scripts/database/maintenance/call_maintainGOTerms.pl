#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance;
use GOTools;

my $dbh = returnProductionDBHNOCommit();
maintainGOTerms($dbh);
$dbh->commit();
