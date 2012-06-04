#! /usr/bin/perl
use strict;
use warnings;
use entrezTools;

if(scalar(@ARGV) != 1){
    print "Usage: call_returnRefFromPMID.pl PMID\n";
    exit;
}

my ($PMID) = @ARGV;;



my $ref=returnRefFromPMID($PMID);
print "PMID: $PMID\n";
foreach my $key (keys %$ref){
	
    print "key: $key\t$ref->{$key}\n";
}

