#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;


if(scalar(@ARGV) != 1){

    print "USAGE: perl addHuman.pl GENE_LIST\n";
}
my $geneList = $ARGV[0];
my $output = $geneList.".human";

open(GENES, $geneList) || die "Cannot open geneList file $geneList for reading\n";
if(!-e $output){
    `touch $output`;
}
open(OUT, ">$output") || die "Can't open output $output file for writing\n";
my $line;
while(defined($line = <GENES>)){
    chomp $line;
    my $uniprot = $line;

    $uniprot .= "_HUMAN";
    print OUT "$uniprot\n";

    
}
