#! /usr/bin/perl
use strict;
use warnings;
use HMMTools;

if(scalar(@ARGV) != 2){
    print "USAGE: call_makeSpeciesFile.pl PFAM_FILE OUTPUT_FILE\n";
    exit;
}

my ($pfamFile, $outputFile) = @ARGV;

makeSpeciesFile($pfamFile, $outputFile);
