#! /usr/bin/perl
use strict;
use warnings;
use compendia;

#retrieve the file, parse the file, combine the parsed files into a single dataset for load, erase the intermediaries.

my $query = "CARBOHYD homo sapiens";
my $uniprotFile = "../data/Uniprot_CARBOHYD.txt";
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
my $dateStr= sprintf("%02d_%02d_%4d", $mon+1, $mday, $year+1900); 
my $outputFile = "../data/Uniprot_CARBOHYD".$dateStr.".txt";

retrieveUniprotFile($query, $uniprotFile);
my $STRICT = 1;
my $numRecordsPerFile = 50000;
my $files = parseModResFile($uniprotFile, $outputFile, $STRICT, $numRecordsPerFile);
print "Files produced @$files\n";
