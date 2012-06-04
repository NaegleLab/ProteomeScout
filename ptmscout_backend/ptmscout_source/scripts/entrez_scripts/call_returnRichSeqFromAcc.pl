#! /usr/bin/perl 
use strict;
use warnings;
use entrezTools;

my $gi;
my $richSeq;

$gi = 'gi|29725609';

$richSeq = returnRichSeqFromAcc($gi);


$gi = 'gi|345600000000000';
$richSeq = returnRichSeqFromAcc($gi);
print "Case 2: Known failure for missing accession\n";
if($richSeq == -1){
    print "PASSED\n";
}
else{
    print "FAILED\n";
}
