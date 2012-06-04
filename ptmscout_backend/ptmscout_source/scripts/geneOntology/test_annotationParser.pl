#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;
use GO::Parser;
use GO::OntologyProvider::OboParser;


# my $annotationFile = "/data/knaegle/data/GO/gene_association.goa_human";
# my $ontologyFile = "/data/knaegle/data/GO/gene_ontology_edit.obo"; #current obo file

# # Try SLIM as well
# $annotationFile = "/home/knaegle/SVN/knaegle/scripts/geneOntology/gene_association.goslim.goa_human";
# #$ontologyFile = "/data/knaegle/data/GO/goslim_generic.obo";

# my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $annotationFile);

# my $ontologyXML = '/data/knaegle/data/GO/go_daily-termdb.obo-xml';

# my $parser = new GO::Parser({handler=>'obj'});
# # my $parser = new GO::Parser({handler=>'obj', use_cache=>1});
# $parser->handler->file("output.xml");
#  $parser->parse($ontologyXML);
#  my $ontology = $parser->handler->graph;

# #my $ontology = GO::OntologyProvider::OboParser->new(ontologyFile => $ontologyFile);

###---------------------------------------------------##

# I want to test annotation parser for species specific -- does this handle uniprotKB and gene names equally
my $annotationFile = "/data/knaegle/data/GO/gene_association.rgd"; # rat
my $annotationParser = GO::AnnotationProvider::AnnotationParser->new(annotationFile => $annotationFile);

my $gene = 'Tdp1';
my $aspect = 'F';
 my $associations = $annotationParser->goIdsByName(name=>$gene, aspect=>$aspect);
print "For Gene: $gene found @$associations\n";
