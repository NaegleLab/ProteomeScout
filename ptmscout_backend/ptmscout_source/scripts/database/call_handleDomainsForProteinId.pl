#! /usr/bin/perl
use strict;
use warnings;
use DBTools::dbEntry;
use DBTools::dbIO;
use HMMTools;

if (scalar(@ARGV) < 1){
    print "Usage: call_handleDomainsForProteinId.pl PROTEIN_ID\n";
    print "INTO TEST DATABASE\n";
    exit;

}

my $dbh = returnTestDBH();
my $proteinId = $ARGV[0];

my $HMM_pvalueCutoff = 0.00001;

my $accArr = returnAccValuesByProteinId($dbh, $proteinId);
my $acc = $accArr->[0];
my ($errorCode, $domainIds) = handleDomainsForProteinId($dbh, $proteinId, $acc, $HMM_pvalueCutoff);
if(!$errorCode){

    print "PASSED\n";
}
else{

    print "FAILED\n";
    print "Error code: $errorCode and domainIDs ".@$domainIds."\n";
}
