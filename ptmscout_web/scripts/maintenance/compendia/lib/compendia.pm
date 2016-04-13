use strict;
use warnings;
use Bio::SeqIO;
use entrezTools;
use commonTools;

# retrieveUniprotFile($query, $outputFile)
# Automatically retrieve the flat text Uniprot file based on the search term
# Inputs: $query - query, for example MOD_RES or PHOSPHORYLATION LARGE SCALE ANALYSIS AT
#         $outputFile - location of file 
# Kristen Naegle
# September 22, 2009
sub retrieveUniprotFile($$){
    my ($query, $outputFile) = @_;
    $query =~ s/ /+/g;
    my $file = "index.html?query=[$query]";
    `rm $file`;
    my $queryString = "http://www.uniprot.org/uniprot/?query=[$query]";
    my $postStr = "--post-data 'format=txt'";
#    &sort=score\&format=txt";
    my $str = $postStr." ".$queryString;
    print "DEBUG: address: $str";
    `wget $str`;
    `mv $file $outputFile`;
}

# parseModResFile($inputFile, $outputFile, $STRICT, $numRecordsPerFile);
# Parses uniprot search for MOD_RES, keeping all modifications and creating output files that are numbered so we can restrict the number of records per file.
# Inputs: $inputFile - txt format search result
#         $outputFile - destination for parsed output
#         $STRICT - boolean, if 1 then only consider MOD_RES known by strict, not by similarity, partial, probably or potential. 
#         $numRecordsPerFile - number of recrods to allow per file
# Kristen Naegle
# Sept. 22, 2009 
# Modified Feb 4, 2013
# Modified Aug 2, 2014 -- Check for probability is over restricitive, catching some probabilities for other terms in the value string
sub parseModResFile($$$$){
    my ($inputFile, $outputFile, $STRICT, $numRecordsPerFile) = @_;
    
    my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');
    my $numberResidues = 15;
    #my $speciesAllowed = returnHashAllowedSpecies();
    my $countLines = 0;
    my $numFile = 1;
    my @files; 
    my $outputFileNew = $outputFile."_$numFile";
    
    open(OUT, ">$outputFileNew") || die "Can't open $outputFileNew for writing\n";
    print OUT "acc\tpep\tMOD_TYPE\n";
    push @files, $outputFileNew;	    

    while(my $seq = $in->next_seq()){
        if($countLines > $numRecordsPerFile){ #check to see if we need to start a new file
            close(OUT);
            $countLines = 0;
            $numFile += 1;
            $outputFileNew = $outputFile."_$numFile";
            open(OUT, ">$outputFileNew") || die "Can't open $outputFile for writing\n";
            print OUT "acc\tpep\tMOD_TYPE\n";
        }
	
        # print $seq."\n";
        #get accession and sequence
        my ($errorCode, $sequence, $species, $gene, $geneSyn, $name, $primaryAcc) = getProteinFromRichSeq($seq);
        my @siteArr;
        for my $feat_object ($seq->get_SeqFeatures){
            if($feat_object->primary_tag eq "MOD_RES"){
            for my $tag ($feat_object->get_all_tags){
                for my $value ($feat_object->get_tag_values($tag)){
                    if($value =~ m/;/){
                        my @value = split(';', $value);
                        $value = $value[0];
                    }
                    if($STRICT){
                        if($value =~ m/(similarity|partial|probable|potential)/i){
                            print "DEBUG: Skipping $value\n";
                            next; #skip those records that are not 
                        }
                    }
                my $pos = $feat_object->location->start;
                my $site = substr($sequence,$pos-1,1);
                my $shortSeq = returnAlignedSequence($numberResidues, $sequence, $site, $pos);
                #print "DEBUG: $shortSeq\n";
                #my $modSeq = $modSeqStart.lc($site).substr($sequence,$pos+2);
                my $modSeq = $seq;
                #return the sequence with the position lower cased.
                my $code = $value; #get the name 
                my $name = $value;
                my $anno = ' ';
                if($value =~ m/./){
                    my @code = split(/\./, $value);
                   $name = $code[0];
                   if(scalar(@code)>1){
                  $anno = $code[1]; 
              }

                }
                my @site = ($shortSeq, $name, $anno);
                push @siteArr, \@site;
		    } #end for my $value
		}
		
	    }
	}
    
	foreach my $site (@siteArr){
        #print "Site info: @$site->[1]\n";
	    my ($modSeq, $name, $anno) = @$site;
	    print OUT "$primaryAcc\t$modSeq\t$name\t$anno\n";
	    $countLines += 1;
	}
	# } #end speices
    
    }
    #foreach my $s (keys %pHash){
#	print OUT "$primaryAcc\t$sequence\t$s\t$pHash{$s}\n";
#    }
    
    close(OUT);
    return \@files;
    
}

# parseCarbohydFile($inputFile, $outputFile, $STRICT, $numRecordsPerFile);
# Parses uniprot search for CARBOHYD, keeping all modifications and creating output files that are numbered so we can restrict the number of records per file.
# Inputs: $inputFile - txt format search result
#         $outputFile - destination for parsed output
#         $STRICT - boolean, if 1 then only consider MOD_RES known by strict, not by similarity, partial, probably or potential. 
#         $numRecordsPerFile - number of recrods to allow per file
# Kristen Naegle
# November 13, 2015 
sub parseCarbohydFile($$$$){
    my ($inputFile, $outputFile, $STRICT, $numRecordsPerFile) = @_;
    
    my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');
    my $numberResidues = 15;
    #my $speciesAllowed = returnHashAllowedSpecies();
    my $countLines = 0;
    my $numFile = 1;
    my @files; 
    my $outputFileNew = $outputFile."_$numFile";
    
    open(OUT, ">$outputFileNew") || die "Can't open $outputFileNew for writing\n";
    print OUT "acc\tpep\tMOD_TYPE\n";
    push @files, $outputFileNew;	    

    while(my $seq = $in->next_seq()){
        if($countLines > $numRecordsPerFile){ #check to see if we need to start a new file
            close(OUT);
            $countLines = 0;
            $numFile += 1;
            $outputFileNew = $outputFile."_$numFile";
            open(OUT, ">$outputFileNew") || die "Can't open $outputFile for writing\n";
            print OUT "acc\tpep\tMOD_TYPE\n";
        }
	
        # print $seq."\n";
        #get accession and sequence
        my ($errorCode, $sequence, $species, $gene, $geneSyn, $name, $primaryAcc) = getProteinFromRichSeq($seq);
        my @siteArr;
        for my $feat_object ($seq->get_SeqFeatures){
            if($feat_object->primary_tag eq "CARBOHYD"){
            for my $tag ($feat_object->get_all_tags){
                for my $value ($feat_object->get_tag_values($tag)){
                    if($value =~ m/;/){
                        my @value = split(';', $value);
                        $value = $value[0];
                    }
                    if($STRICT){
                        if($value =~ m/(similarity|partial|probable|potential)/i){
                            print "DEBUG: Skipping $value\n";
                            next; #skip those records that are not 
                        }
                    }
                my $pos = $feat_object->location->start;
                my $site = substr($sequence,$pos-1,1);
                my $shortSeq = returnAlignedSequence($numberResidues, $sequence, $site, $pos);
                print "DEBUG: $shortSeq\n";
                #my $modSeq = $modSeqStart.lc($site).substr($sequence,$pos+2);
                my $modSeq = $seq;
                #return the sequence with the position lower cased.
                my $code = $value; #get the name 
                my $name = $value;
                my $anno = ' ';
                if($value =~ m/./){
                    my @code = split(/\./, $value);
                   $name = $code[0];
                   if(scalar(@code)>1){
                  $anno = $code[1]; 
              }

                }
                my @site = ($shortSeq, $name, $anno);
                push @siteArr, \@site;
		    } #end for my $value
		}
		
	    }
	}
    
	foreach my $site (@siteArr){
        #print "Site info: @$site->[1]\n";
	    my ($modSeq, $name, $anno) = @$site;
	    print OUT "$primaryAcc\t$modSeq\t$name\t$anno\n";
	    $countLines += 1;
	}
	# } #end speices
    
    }
    #foreach my $s (keys %pHash){
#	print OUT "$primaryAcc\t$sequence\t$s\t$pHash{$s}\n";
#    }
    
    close(OUT);
    return \@files;
    
}




# $alignedShortenedSequence = returnAlignedSequence($numAA, $sequence, $code, $alignNum);
# Returns the shortened sequence (2*$numAA + 1 in length) from the full sequence aligned on the alignment number and with code (S, T, or Y) 
# Inputs: $numAA the number of amino acids to the n-terminal and c-terminal side
#         $sequence the full AA sequence 
#         $code can be S, T, or Y represents the phosphorylated residue
#         $alignNum is the position in the full sequence (ones-based) of the $code 
# Test file: compareDataTopELM.pl
# Kristen Naegle, 4/12/07
sub returnAlignedSequence($$$$){
    my $numAA = shift;
    my $strSequence = shift;
    my $code = shift;
    my $alignNum = shift;

    my @sequence = split(//, $strSequence);
    my $pos = $alignNum - 1; #shift to 0-based counting
    
    #if the position and the code are aligned then return the substring
    if($sequence[$pos] eq $code){
	my $posStart = $pos - $numAA;
	my $numStart = $numAA;
	my $numNPad = 0;
	my $numCPad = 0;
	#handle N-terminus
	if ($posStart < 0) {
	    $numNPad = $numAA - $pos;
	    $numStart = -$posStart;
	    $posStart = 0;
	}
	my $posEnd = $pos + $numAA;
	#handle C-terminus
	if ($posEnd > $#sequence){
	    $numCPad = $posEnd - $#sequence;
	    #print("Need to pad $numCPad where $posEnd is desired end and $#sequence is C terminmal\n"); 
	    $posEnd = $#sequence;

	}
	$sequence[$pos] = lc($sequence[$pos]); #if phosphorylated change to a lowercase 
	my @subSeq = @sequence[$posStart..$posEnd];
	my $subStrSeq = join("", @subSeq);
	if($numNPad){
	    my $i;
	    for($i=0; $i < $numNPad; $i++){
		$subStrSeq = " ".$subStrSeq;
	    }
	}
	if($numCPad){
	    my $i;
	    for($i=0; $i < $numCPad; $i++){
		$subStrSeq = $subStrSeq." ";
	    }
	}
	return $subStrSeq;
    }
    # else if they don't align return zero
    else {
	print "In sub returnAlignedSequences $sequence[$pos] for position:$pos does not equal $code\n";
	return 0;
    }
}

1;
