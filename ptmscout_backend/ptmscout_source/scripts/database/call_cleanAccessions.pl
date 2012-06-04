#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use DBTools::maintenance;
use errorHandling;
use HMMTools;

my $dbh = returnProductionDBHNOCommit();
flushLog();

cleanAccessions($dbh);

$dbh->commit();
