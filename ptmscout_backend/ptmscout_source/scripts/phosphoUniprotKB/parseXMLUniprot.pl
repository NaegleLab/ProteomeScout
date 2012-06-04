#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;


#my $inputFile = 'uniprot-PHOSPHORYLATION-LARGE-SCALE-ANALYSIS-AT.txt';
my $inputFile = 'smallResultsFlat.txt';

my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');

my $search = "PHOSPHORYLATION \[LARGE SCALE ANALYSIS\] AT";

while(my $seq = $in->next_seq()){
    print $seq."\n";

    # show that I can print the accession number (first one) and scope terms for Phosphorylation at terms:
#    my ($ann) = $seq->annotation->get_Annotations('gene_name');
#    my @names = $ann->findval('Name');
#    print "Names: @names\n";
#     foreach my $key (keys %$seq){
# 	print "\tKey: $key\n";
# 	print "\t$seq->{$key}\n";
#     }

    print "Annotation Collection\n";
    my $ac = $seq->annotation();
    foreach my $k ($ac->get_all_annotation_keys()){
	my @values = $ac->get_Annotations($k);
	foreach my $value (@values){
	    #print "Annotation ", $k, "stringfiled value ",$value->as_text,"\n";
	    if($value->tagname eq "reference"){
		my $rp = $value->rp();
		if($rp =~ m/PHOSPHORYLATION/ and $rp =~ m/LARGE SCALE ANALYSIS/){
#		my $hash_ref = $value->hash_tree();
#		for my $key (keys %$hash_ref){
#		    print $key,": ",$hash_ref->{$key},"\n";
#		}
		    print "RP: ",$value->rp(),"\n";
	    
		}
	    }

	}
    }

    print "Features\n";
    for my $feat_object ($seq->get_SeqFeatures) {          
	print "primary tag: ", $feat_object->primary_tag, "\n";          
	if($feat_object->primary_tag eq "MOD_RES"){
	for my $tag ($feat_object->get_all_tags) {             
	    print "  tag: ", $tag, "\n";         

		for my $value ($feat_object->get_tag_values($tag)) {                
		    print "    value: ", $value, "\n";     
		    print "         start: ",$feat_object->location->start,"\n";
		}          
	    }
	}       
    }

    my @ModRes;
    for my $feat_object ($seq->get_SeqFeatures) {    
	push @ModRes, $feat_object->get_tag_values("Phosphoserine") if ($feat_object->has_tag("Phosphoserine"));     
    }
    print "FOUND ".scalar(@ModRes)." modified serine residues\n";

}
