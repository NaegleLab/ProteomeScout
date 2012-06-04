#! /usr/bin/perl
use strict;
use warnings;
use GO::AnnotationProvider::AnnotationParser;
use DBTools::dbIO;
use DBTools::queryTools;
use DBTools::insertFunctions;
use GOTools;
use DBTools::dbEntry;
use errorHandling;
use GO::Parser;


if(scalar(@ARGV) != 2){

    print "USAGE: perl call_handleGoTermsForProteinIds_viaExpID.pl EXPERIMENT_ID SLIM\n";
    exit;
}

flushLog();
#my $geneList = $ARGV[0];
my $expId = $ARGV[0];
my $SLIM = $ARGV[1];
#my $SLIM = 1;

my $dbh = returnProductionDBHNOCommit();

handleGoTermsForProteinIds_viaExpId($dbh, $expId, $SLIM);


$dbh->commit();
