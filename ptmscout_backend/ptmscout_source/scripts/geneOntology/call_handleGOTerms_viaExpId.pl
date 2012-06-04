#! /usr/bin/perl
use strict;
use warnings;
use GOTools;
use DBTools::dbEntry;
use DBTools::dbIO;

if(scalar(@ARGV) != 1){
    print "Usage: call_handleGOTermsForProteinIds_viaExpId.pl EXP_ID\n";
    exit;
}

my($expId) = @ARGV;
my $dbh = returnProductionDBHNOCommit();

my ($ontologyFile);

$ontologyFile = "/home/knaegle/SVN/knaegle/scripts/GO_DATA/gene_ontology.1_2.obo";


# Make sure that proteinIds with other gene_synonyms is found as well

my $GOHash = handleGOTermsForProteinIds_viaExpId($dbh,$expId, $ontologyFile);
my @proteins = keys(%$GOHash);
print "@proteins \n ".scalar(@proteins)." successfully had GO terms inserted\n";
$dbh->commit();
