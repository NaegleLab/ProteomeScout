#! /usr/bin/perl
use strict;
use warnings;
use dataIO;

if(@ARGV != 2){
    print "USAGE: call_printRefFromPMID.pl PMID OUTPUT_FILE\n";
    exit;
    
}

my($PMID, $outputFile) = @ARGV;

printRefFromPMID($PMID, $outputFile);

