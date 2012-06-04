#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::maintenance; 
use DBTools::resultParsing;

# This will update all pfam_site in phosphopep table for the database

my $dbh = returnProductionDBHNOCommit();


# get all phosphopeps
my $sth = $dbh->prepare('SELECT id from phosphopep');
$sth->execute();
my $phosphopeps = returnArrayOfResultsOnCol($sth,0);

updatePfamSitesForPhosphopeps($dbh, $phosphopeps);


$dbh->commit;
