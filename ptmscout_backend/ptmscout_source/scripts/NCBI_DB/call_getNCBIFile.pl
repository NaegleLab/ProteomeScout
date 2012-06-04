#! /usr/bin/perl
use strict;
use warnings;
use ncbiDBTools;


if(scalar(@ARGV) != 1){

    print "Usage: call_getNCBIFile.pl FILE_NUMBER\n";
    exit;
}
my $number = $ARGV[0];
getNCBIFile($number);



