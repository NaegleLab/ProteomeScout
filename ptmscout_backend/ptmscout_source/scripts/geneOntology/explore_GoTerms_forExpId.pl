#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;
use DBTools::dbIO;
use DBTools::queryTools;
use DBTools::insertFunctions;
use GOTools;


if(scalar(@ARGV) != 3){

    print "USAGE: perl explore_GoTerms.pl EXPERIMENT_ID <C|F|P> BOOL_SLIM\n";
    exit;
}
#my $geneList = $ARGV[0];
my $expId = $ARGV[0];
my $aspect = $ARGV[1];
my $SLIM = $ARGV[2];

my $dbh = returnProductionDBHNOCommit();

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

#open(GENES, $geneList) || die "Cannot open geneList file $geneList for reading\n";
my $accHash = returnGeneListForExperiment($dbh, $expId);

#my $uniprot = 'P00533';
#my $line;
#while(defined($line = <GENES>)){
#my $aspect = 'P';
foreach my $protein (keys %$accHash){
    print "\n";
    my $uniprot = $accHash->{$protein};
    print "Protein ID: $protein\t Acc: $uniprot\n";
   
 # if($uniprot eq '-1'){
# 	print "ERROR: No accession\n";
# 	next;
#     }
#     print "ID: $uniprot\n";
#     if($annotationParser->nameIsAnnotated(name=>$uniprot)){
# 	print "YES $uniprot\n";
#     }
#     else{
# 	print "NO $uniprot\n";
# 	my $name;
# 	($uniprot, $name) = returnStandardNameForAlias($annotationFile, $uniprot, '9606');
# 	if($uniprot eq '-1'){
# 	    print "ERROR: Could not find alias $uniprot\n";
# 	    next;
# 	}
# 	else{
# 	    # here can insert new swissprot identity for that protein
# 	}
#     }
#     my $associations = $annotationParser->goIdsByName(name=>$uniprot, aspect=>$aspect);
#     if(not defined($associations)){
# 	print "NOT DEFINED...Looking for it in standard name\n";
# 	#$uniprot .= "_HUMAN";
# 	if($annotationParser->nameIsAmbiguous($uniprot)){
# 	    print "ERROR: Ambiguous!! \n";
# 	    next;
# 	}
# 	my $associations = $annotationParser->goIdsByName(name=>$uniprot, aspect=>$aspect);
# 	if(not defined($associations)){
# 	    print "STILL Not found by standard name..\n";
# 	    my $dbId = $annotationParser->databaseIdByName(name=>$uniprot, aspect=>$aspect);
# 	    print "CHECKED for dbID: $dbId\n";

# 	}
#     }
    my ($error, $associations) = returnGOAssociations($annotationParser, $annotationFile, $uniprot, $aspect, $taxonomy);
    if($error){
	print "ERROR -- unable to find association for gene: $uniprot\n";
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
