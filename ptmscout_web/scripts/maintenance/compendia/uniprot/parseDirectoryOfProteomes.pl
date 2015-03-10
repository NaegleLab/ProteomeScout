#! /usr/bin/perl
use lib "../lib";
use strict;
use warnings;
use compendia;

#retrieve the file, parse the file, combine the parsed files into a single dataset for load, erase the intermediaries.

if(scalar(@ARGV) != 2){
	print "Usage: parseModResFile.pl MOD_RES_DIRECTORY NUM_RECSPERFILE\n";
	exit;
}
my ($uniprotDir, $numRecordsPerFile) = @ARGV;
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdist) = localtime(time);
my $dateStr= sprintf("%02d_%02d_%4d", $mon+1, $mday, $year+1900); 
#my $outputFile = "../data/Uniprot_MODRES".$dateStr.".txt";
my $outputDir = $uniprotDir.'output/'; 
my $outputFile = $outputDir.'Uniprot_MODRES'.$dateStr.".txt";

my $STRICT = 1;
#create one master file to send into parser and create from that
my $masterFile = $outputDir.'Uniprot_MODRES_master.txt'; 
open(OUT, ">$masterFile") || die "Can't open $masterFile for writing\n";
foreach my $fp (glob("$uniprotDir/*.txt")){
    open(FILE, $fp) || ((warn "Can't open file $fp\n"), next FILE);
    while (<FILE>){
        print OUT;
    }
    close(FILE);

}
my $files = parseModResFile($masterFile, $outputFile, $STRICT, $numRecordsPerFile);
print "Files produced @$files\n";
