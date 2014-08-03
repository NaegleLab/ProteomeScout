#! /usr/bin/perl
use strict;
use warnings;
use compendia;

#retrieve the file, parse the file, combine the parsed files into a single dataset for load, erase the intermediaries.

if(scalar(@ARGV) != 2){
	print "Usage: parseModResFile.pl MOD_RES_FILE NUM_RECSPERFILE\n";
	exit;
}
my ($uniprotFile, $numRecordsPerFile) = @ARGV;
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
my $dateStr= sprintf("%02d_%02d_%4d", $mon+1, $mday, $year+1900); 
my $outputFile = "../data/Uniprot_MODRES".$dateStr.".txt";

#retrieveUniprotFile($query, $uniprotFile);
my $STRICT = 1;
my $files = parseModResFile($uniprotFile, $outputFile, $STRICT, $numRecordsPerFile);
print "Files produced @$files\n";
