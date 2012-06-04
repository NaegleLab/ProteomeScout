#! /usr/bin/perl
use strict;
use warnings;
use DBI;
use DBTools::dbEntry;
use DBTools::dbIO;
use DBTools::insertFunctions;


print "MAKE SURE TO DELETE CURRENT PELM ANNOTATIONS before calling this script\n";

my $dbh = returnProductionDBHNOCommit();

flushLog();
my $pelmFile = "/data/knaegle/data/phosphoELM/phosphoELM_current.txt";
my($RELEASE, $file) = returnLinkRelease($pelmFile);
$RELEASE =~ s/\.txt//;

#my $release_date = "2008";
my $predIds;


$predIds = handlePELMKinaseAnnotations($dbh, $pelmFile);
print "Inserted ".scalar(@$predIds)." into database from PELM\n";

$dbh->commit();
