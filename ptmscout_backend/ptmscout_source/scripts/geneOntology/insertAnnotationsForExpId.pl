#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;
use DBTools::dbIO;
use DBTools::queryTools;
use DBTools::insertFunctions;
use GOTools;


if(scalar(@ARGV) != 2){

    print "USAGE: perl explore_GoTerms.pl EXPERIMENT_ID <C|F|P> \n";
    exit;
}
#my $geneList = $ARGV[0];
my $expId = $ARGV[0];
my $aspect = $ARGV[1];
#my $SLIM = $ARGV[2];
my $SLIM = 1;

my $dbh = returnTestDBHNOCommit();

my $taxonomy = '9606';
my $annotationParser;
my $annotationFile;
if($SLIM){
    $annotationFile = "/home/knaegle/SVN/knaegle/scripts/geneOntology/gene_association.goslim.goa_human";
}
else{
    $annotationFile = "/data/knaegle/data/GO/gene_association.goa_human";

}
    $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $annotationFile);

my $accHash = returnGeneListForExperiment($dbh, $expId);

foreach my $protein (keys %$accHash){
    print "\n";
    my $uniprot = $accHash->{$protein};
    print "Protein ID: $protein\t Acc: $uniprot\n";
    my ($error, $associations) = returnGOAssociations($annotationParser, $annotationFile, $uniprot, $aspect, $taxonomy);
    if(!$error){
	#1. check to see of GO is already in database, if not insert both tables
	#2. if so then check to see if protein is already linked to ontology .. if not then add this table
	
    }
#    print "GO Associations for gene: \n";
    foreach my $a (@$associations){
	print "\t$a\n";
    }

#     print "Database ID for gene: ", $annotationParser->databaseIdByName($uniprot), "\n";

#     print "Database name: ", $annotationParser->databaseName(), "\n";    
#     print "Standard name for gene: ", $annotationParser->standardNameByName($uniprot), "\n";

}


$dbh->rollback;
