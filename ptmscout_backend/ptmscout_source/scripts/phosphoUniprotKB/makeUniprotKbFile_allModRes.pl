#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use entrezTools;

# This takes in a flat text format of a search or the entire Uniprot database and outputs sequence, site:pos, site:code and primary accession number.  This will look in all MOD_RES features for Phospho-serine/threonine/and tyrosine features.  

if(scalar(@ARGV) != 3){
    print "USAGE: makeUniprotKbFile_allModRes.pl INPUT_FILE OUTPUT_FILE STRICT\n";
    print "INPUT_FILE is flat file of uniprot form, OUTPUT_FILE is target text and STRICT is Boolean meaning to 1- restrict search to non-partial and not by homology sites, 0, include all\n";
    exit;
}

my ($inputFile, $outputFile, $STRICT) = @ARGV;
open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');

my %transhash;
$transhash{'Phosphoserine'} = 'S';
$transhash{'Phosphothreonine'} = 'T';
$transhash{'Phosphotyrosine'} = 'Y';

print OUT "acc:sp\tsequence:full\tsite:pos\tsite:code\n";
while(my $seq = $in->next_seq()){
   # print $seq."\n";
    #get accession and sequence
    my ($errorCode, $sequence, $species, $gene, $geneSyn, $name, $primaryAcc) = getProteinFromRichSeq($seq);
    my %siteHash;
    
    # push all the phosphorylations onto a hash then print
    for my $feat_object ($seq->get_SeqFeatures){

	if($feat_object->primary_tag eq "MOD_RES"){
	    for my $tag ($feat_object->get_all_tags){
	#	print "tag: ", $tag, "\n";
		for my $value ($feat_object->get_tag_values($tag)){
		    my $match = '(Phosphoserine|Phosphotyrosine|Phosphothreonine)';
		    my $INSERT = 1;
		    if($value =~ m/$match/){
			my $code = $transhash{$1};
			my $site = $feat_object->location->start;

			if($STRICT){
			    if($value =~ m/(similarity|partial|probably|potential)/i){
				$INSERT = 0;
			    }
			}
			if(not defined $siteHash{$site} and $INSERT){
			    $siteHash{$site}=$code;
			}
		    }
		   # print "    value: ", $value, "\n";
		   # print "         start: ", $feat_object->location->start, "\n";
		}
	    }

	}
    }
    foreach my $s (keys %siteHash){
	print OUT "$primaryAcc\t$sequence\t$s\t$siteHash{$s}\n";
    }
	
}
    #foreach my $s (keys %pHash){
#	print OUT "$primaryAcc\t$sequence\t$s\t$pHash{$s}\n";
#    }

close(OUT);
