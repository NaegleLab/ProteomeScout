#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use entrezTools;

# This takes in a flat text format of a search or the entire Uniprot database and outputs sequence, site:pos, site:code and primary accession number.  Only for those things indicated by PHOSPHORYLATION [LARGE SCALE ANALYSIS] AT

if(scalar(@ARGV) != 2){
    print "USAGE: makeUniprotKbFile_largeScaleOnly.pl INPUT_FILE OUTPUT_FILE\n";
    exit;
}

my ($inputFile, $outputFile) = @ARGV;
#my $inputFile = 'uniprot-PHOSPHORYLATION-LARGE-SCALE-ANALYSIS-AT.txt';
#my $inputFile = 'smallResultsFlat.txt';
#my $outputFile = 'smallResults_out.txt';
open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');

my $search = "PHOSPHORYLATION [LARGE SCALE ANALYSIS] AT";
print OUT "acc:sp\tsequence:full\tsite:pos\tsite:code\n";
while(my $seq = $in->next_seq()){
   # print $seq."\n";
    #get accession and sequence
    my ($errorCode, $sequence, $species, $gene, $geneSynonyms, $name, $primaryAcc) = getProteinFromRichSeq($seq);
    my %siteHash;
    
    # push all the phosphorylations onto a hash then print
    my %pHash;
    my $ac = $seq->annotation();
    foreach my $k ($ac->get_all_annotation_keys()){
	my @values = $ac->get_Annotations($k);
	foreach my $value (@values){
	    if($value->tagname eq "reference"){
		my $rp = $value->rp();
		if($rp =~ m/PHOSPHORYLATION/ and $rp =~ m/LARGE SCALE ANALYSIS/){
		    my ($sites, $codes) = parseSiteCodesFromRP($rp);
		    #print "Found ".scalar(@$sites)."\n";
		    #print "RP: ",$value->rp(),"\n";
		    for(my $i=0; $i<scalar(@$sites); $i++){
			my $site = $sites->[$i];
			my $code = $codes->[$i];
			#print "Site: $site\t Code: $code\n";
			if(not defined($pHash{$site})){
			    $pHash{$site} = $code;
			}
		    }
	    
		}
	    }

	} #end of every reference so now print sites

	
    }
    foreach my $s (keys %pHash){
	print OUT "$primaryAcc\t$sequence\t$s\t$pHash{$s}\n";
    }
}
close(OUT);
