#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::deleteFunctions;


my $dbh = returnProductionDBHNOCommit();

clearAllNonExpressionTables($dbh);

$dbh->commit();
