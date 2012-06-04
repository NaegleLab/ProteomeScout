#! /usr/bin/perl
use strict;
use warnings;
use expressionTools;
use DBTools::dbIO;

my $dbh = returnProductionDBHNOCommit();
my ($numHandled, $numFailed);

($numHandled, $numFailed) = maintainExpressionLinkages($dbh);
print "Handled: $numHandled\tFAILED: $numFailed\n";


$dbh->commit();
