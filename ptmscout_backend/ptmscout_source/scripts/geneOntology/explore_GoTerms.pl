#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;


if(scalar(@ARGV) != 2){

    print "USAGE: perl explore_GoTerms.pl GENE_LIST BOOL_SLIM\n";
}
my $geneList = $ARGV[0];
my $SLIM = $ARGV[1];

my $annotationParser;
my $annotationFile;
if($SLIM){
    $annotationFile = "/home/knaegle/SVN/knaegle/scripts/geneOntology/gene_association.goslim.goa_human";
}
else{
    $annotationFile = "/data/knaegle/data/GO/gene_association.goa_human";

}
    $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $annotationFile);

open(GENES, $geneList) || die "Cannot open geneList file $geneList for reading\n";

#my $uniprot = 'P00533';
my $line;
while(defined($line = <GENES>)){
    chomp $line;
    my $uniprot = $line;

   # $uniprot .= "_HUMAN";
    print "ID: $uniprot\n";
#my $uniprot = 'EGFR_HUMAN';
    my $associations = $annotationParser->goIdsByName(name=>$uniprot, aspect=>'P');
    print "GO Associations for gene: \n";
    foreach my $a (@$associations){
	print "\t$a\n";
    }

    print "Database ID for gene: ", $annotationParser->databaseIdByName($uniprot), "\n";

    print "Database name: ", $annotationParser->databaseName(), "\n";    
    print "Standard name for gene: ", $annotationParser->standardNameByName($uniprot), "\n";
    my $i;
    my @geneNames = $annotationParser->allStandardNames();
    foreach $i (0..10) {    
#	print "$geneNames[$i]\n";    
    }
}
