#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use errorHandling;
use HMMTools;

if(scalar(@ARGV) != 1){
    print "Usage: call_loadProteinsAndDomainsInFile.pl DATA_FILE \n";
    exit;
}

my $dataFile = $ARGV[0];


my $dbh = returnProductionDBHNOCommit();
flushLog();

my ($expId, $accepted, $rejected);

#($expId, $accepted, $rejected) = loadDataFile($dbh, $dataFile,$name, $author, $description);
my $hash = loadProteinAndDomainsInFile($dbh, $dataFile);

print "Accession\tproteinId\tiProteinFlag\tiAccFlag\n";
foreach my $k (keys %$hash){
    my @array = @{$hash->{$k}};
    print "Accession:\t$k\t@array\n";
}

# Time stamp log.
my $newLog = dateStampLogFile();
print "LOG AT: $newLog\n";   

$dbh->commit();
