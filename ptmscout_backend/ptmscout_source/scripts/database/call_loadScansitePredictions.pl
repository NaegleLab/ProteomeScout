#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::insertFunctions;
use DBTools::dbIO;
use DBTools::dbEntry;

my $dbh = returnProductionDBHWithCommit();

loadScansitePredictions($dbh);

$dbh->commit();
#$dbh->rollback();
