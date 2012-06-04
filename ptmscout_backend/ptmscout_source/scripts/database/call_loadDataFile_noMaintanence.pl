#! /usr/bin/perl 
use strict;
use warnings;
use DBTools::dbIO;
use DBTools::dbEntry;
use errorHandling;
use HMMTools;

if(scalar(@ARGV) != 13){
    print "Usage: call_loadDataFile.pl DATA_FILE NAME AUTHOR DESCRIPTION CONTACT PMID URL PUBLISHED AMBGUITY EXPORT EXPID_LINK SUBMISSION_EMAIL PRIMARY_MOD(S)\n";
    exit;
}

my $dataFile = $ARGV[0];
my $name = $ARGV[1];
my $author = $ARGV[2];
my $description = $ARGV[3];
my $contact = $ARGV[4];
my $PMID = $ARGV[5];
my $URL = $ARGV[6];
my $published = $ARGV[7];
my $AMBIGUITY = $ARGV[8];
my $EXPORT = $ARGV[9];
my $expId_link = $ARGV[10];
my $submissionEmail = $ARGV[11];
my $primaryMods = $ARGV[12];

my $dbh = returnProductionDBHNOCommit();
#xsmy $dbh = returnWebDevDBHNOCommit();
#my $dbh = returnTestDBHNOCommit();
flushLog();

my ($expId, $accepted, $rejected);

($expId, $accepted, $rejected) = loadDataFile($dbh, $dataFile,$name, $author, $description, $contact, $PMID, $URL, $published, $AMBIGUITY, $EXPORT, $expId_link, $submissionEmail, $primaryMods);
print "Experiment ID: $expId\n";
print "Accepted: $accepted Rejected:$rejected\n";


