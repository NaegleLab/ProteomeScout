#! /usr/bin/perl 
use strict;
use warnings;
use dataIO;

if(scalar(@ARGV) != 2){
    print "USAGE: addSpacesToPepFile.pl INPUT_FILE OUTPUT_FILE\n";
    exit;
}

my ($inputFile, $outputFile) = @ARGV;

open(IN, $inputFile) || die "Can't open $inputFile for reading\n";
open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";

while(defined(my $line = <IN>)){
    chomp $line;
    my $spaced = returnAlignedSpace($line);
    print OUT $spaced."\n";
    
}

