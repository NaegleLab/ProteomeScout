#! /usr/bin/perl 
use strict;
use warnings;
use GOTools;

if(scalar(@ARGV) != 1){
    print "Usage: call_retrieveAllGOAnnotationFiles.pl REMOVE_IEA\n";
    exit;
}

my($REMOVE_IEA) = @ARGV;

retrieveAllGOAnnotationFiles($REMOVE_IEA);
