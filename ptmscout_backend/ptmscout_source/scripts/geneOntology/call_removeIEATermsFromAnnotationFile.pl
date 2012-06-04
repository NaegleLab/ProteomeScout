#! /usr/bin/perl 
use strict;
use warnings;
use GOTools;

if(scalar(@ARGV) != 2){
    print "USAGE call_removeIEATermsFromAnnotationFile INPUT_FILE OUTPUT_FILE\n";
    exit;
}

my ($inputFile, $outputFile) = @ARGV;

removeIEATermsForAnnotationFile($inputFile, $outputFile);
