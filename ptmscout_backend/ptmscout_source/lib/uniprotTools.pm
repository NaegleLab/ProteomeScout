use strict;
use warnings;
use Bio::SeqIO;
use entrezTools;
use commonTools;



# TO ADD functionality for a new modification, make sure to add it to query, parsing and parseSiteCodesFromRP. 

# retriveUniprotFile($query, $outputFile)
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

# parseLargeScaleFile($inputFile, $outputFile)
# Parse a LARGE SCALE ANALYSIS AT .txt search from uniprot
# Inputs: $modification - the keyword to search for in parsing. (e.g. PHOSPHORYLATION or ACETYLATION).
#         $inputFile - the .txt search report from Uniprot
#         $outputFile - the destination .txt file that contains sequence, site and code information about the modifications
# Kristen Naegle
# Sept. 22, 2009
sub parseLargeScaleFile($$$){
    my ($modification, $inputFile, $outputFile) = @_;
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');

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
		    if($rp =~ m/$modification/ and $rp =~ m/LARGE SCALE ANALYSIS/){

			my ($sites, $codes) = parseSiteCodesFromRP($modification, $rp);
		#	print "DEBUG: met criteria $modification and codes are @$codes\n";
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
    close(OUT)
	
    }

# parseModResFile($modification, $inputFile, $outputFile, $STRICT);
# Parses uniprot search for MOD_RES, keeping all allowed modifications
# Inputs: $inputFile - txt format search result
#         $outputFile - destination for parsed output
#         $STRICT - boolean, if 1 then only consider MOD_RES known by strict, not by similarity, partial, probably or potential. 
# Kristen Naegle
# Sept. 22, 2009 
sub parseModResFile($$$){
    my ($inputFile, $outputFile, $STRICT) = @_;
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    my $in = Bio::SeqIO->new(-file => $inputFile, -format => 'swiss');
    
    my $speciesAllowed = returnHashAllowedSpecies();

    my %transhash;
    $transhash{'Phosphoserine'} = 'S';
    $transhash{'Phosphothreonine'} = 'T';
    $transhash{'Phosphotyrosine'} = 'Y';
    $transhash{'acetyllysine'} = 'K';
    
    my $match = '(Phosphoserine|Phosphotyrosine|Phosphothreonine|acetyllysine)';
    print OUT "acc:sp\tsequence:full\tsite:pos\tsite:code\n";
    while(my $seq = $in->next_seq()){
	# print $seq."\n";
	#get accession and sequence
	my ($errorCode, $sequence, $species, $gene, $geneSyn, $name, $primaryAcc) = getProteinFromRichSeq($seq);
	my %siteHash;
	#
	if($speciesAllowed->{$species}){
	    #print "DEBUG: Species: $species\n";
	    # push all the phosphorylations onto a hash then print
	    for my $feat_object ($seq->get_SeqFeatures){

		if($feat_object->primary_tag eq "MOD_RES"){
		    for my $tag ($feat_object->get_all_tags){
			#	print "tag: ", $tag, "\n";
			for my $value ($feat_object->get_tag_values($tag)){
			    
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
	}
	   foreach my $s (keys %siteHash){
	       print OUT "$primaryAcc\t$sequence\t$s\t$siteHash{$s}\n";
	   }
    }
       
       
    #foreach my $s (keys %pHash){
#	print OUT "$primaryAcc\t$sequence\t$s\t$pHash{$s}\n";
#    }
    
    close(OUT);
    
    
}

# ($phosphoLS, $acetylLS, $modRes, $modRes_STRICT) = returnGlobalUniprotFileLocations();
# Returns the locations of the global files for uniprot. 
# Outputs: $phosphoLS - phosphorylation at large scale analysis 
#          $acetylLS - acetylation Large scale analysis
#          $modRes - MOD_RES search
#          $modRes - MOD_RES search parsed with STRICT
# Kristen Naegle
# Sept. 22, 2009
sub returnGlobalUniprotFileLocations(){
   # my $rootDir = '/data/knaegle/data/uniprot/';
    my $rootDir = "./";
    my $phosphoLS = $rootDir."uniprot-PHOSPHORYLATION-LARGE-SCALE-ANALYSIS-AT.txt";
    my $acetylLS = $rootDir."uniprot-ACETYLATION-LARGE-SCALE-ANALYSIS-AT.txt";
    my $modRes = $rootDir."uniprot-MOD_RES.txt";
    my $modRes_STRICT = $rootDir."uniprot-MOD_RES-STRICT.txt";

    return ($phosphoLS, $acetylLS, $modRes, $modRes_STRICT);

}

# retriveAndParseUniprotFileSearches()
# retrieves, parses accordingly, and then appends the 15-mer fragment to all Uniprot searches (PHOSPHORYLATION, ACETYLATIOn and MOD_RES)
# Kristen Naegle
# September 22, 2009
sub retrieveAndParseUniprotFileSearches(){
    my ($phosphoLS, $acetylLS, $modRes, $modRes_STRICT) = returnGlobalUniprotFileLocations();
    my $date = returnMonthYear();
    my $dataRoot = $globalVars::DATA_PATH;
    my $targetFile;
    my $datedFile;
    #retrieve LARGE SCALE ANALYSIS and then parse them. 
    my $query = "PHOSPHORYLATION LARGE SCALE ANALYSIS AT";
    my $outputFile = 'temp.txt'; # target of initial download
    my $tempFile = 'temp2.txt'; #target of parsing
    retrieveUniprotFile($query, $outputFile);
    parseLargeScaleFile('PHOSPHORYLATION', $outputFile, $tempFile);
    append15merFragment($tempFile, $phosphoLS);
    $targetFile = $dataRoot.'uniprot-PHOSPHORYLATION-LARGE-SCALE';
    $datedFile = $targetFile."-".$date;
    moveAndLink($phosphoLS, $targetFile, $datedFile);

    $query = "ACETYLATION LARGE SCALE ANALYSIS AT";
    retrieveUniprotFile($query, $outputFile);
    parseLargeScaleFile('ACETYLATION', $outputFile, $tempFile);
    print "DEBUG: just parsed into $tempFile\n";
    append15merFragment($tempFile, $acetylLS);
    $targetFile = $dataRoot.'uniprot-ACETYLATION-LARGE-SCALE';
    $datedFile = $targetFile."-".$date;
    moveAndLink($acetylLS, $targetFile, $datedFile);

    $query = "MOD_RES";
    my $STRICT = 0;
    retrieveUniprotFile($query, $outputFile);
    parseModResFile($outputFile, $tempFile, $STRICT);
    append15merFragment($tempFile, $modRes);
    $targetFile = $dataRoot."uniprot-MOD_RES";
    $datedFile = $targetFile."-".$date;
    moveAndLink($modRes, $targetFile, $datedFile);

    $query = "MOD_RES";
    # copy last since missing file
#    `cp $datedFile $tempFile`;
    $STRICT = 1;
    parseModResFile($outputFile, $tempFile, $STRICT);
    append15merFragment($tempFile, $modRes_STRICT);
     $targetFile = $dataRoot."uniprot-MOD_RES_STRICT";
    $datedFile = $targetFile."-".$date;
    moveAndLink($modRes_STRICT, $targetFile, $datedFile);
}

# moveAndLink($outputFile, $targetFile, $datedFile);
# Given a temporary file, move it to it's dated target file and forcibly update symbolic link of non-dated file.
# Inputs: $outputFile - this is the temp file that will be moved to $datedFile
#         $targetFile - name of symbolic link
#         $datedFile - target move of output temp file
# Kristen Naegle
# Nov 25, 2009
sub moveAndLink($$$){
    my($outputFile, $targetFile, $datedFile) = @_;
    `mv $outputFile $datedFile`;
    `ln -sf $datedFile $targetFile`;

}

# Given an input file with site:code, sequence, and site:type, turn this into a file that can be loaded, add 15-mer peptide to it.
# Inputs: $inputFile - .txt input file with above columns
#         $outputFIle - destination file
# Kristen Naegle
# Sept. 22, 2009
sub append15merFragment($$){
    my ($dataFile, $outFile) = @_;
    if(!-e $outFile){`touch $outFile`;}
    open(OUT, ">$outFile") || die "Can't open output file $outFile for writing\n";

    my $seqCol = returnColumnNumber($dataFile, "sequence");
    my $posCol = returnColumnNumber($dataFile, "site:pos");
    my $codeCol = returnColumnNumber($dataFile, "site:code");
    my @colNameArr = ["sequence", "site:pos", "site:code"]; 
    my @colNumArr = [$seqCol, $posCol, $codeCol];
    my $throw = catchColumnError(\@colNameArr, \@colNumArr, "append15merFragment");
    if($throw){
	exit;
    }

    open(DATA, $dataFile) || die "Can't open file $dataFile for reading\n";
    my $line = <DATA>; #skip header
    chomp $line;
    $line .= "\tpep\n";
    print OUT $line;
    my $lineNum = 1;
    while(defined($line = <DATA>)){

	chomp $line;
	my @line = split("\t",$line); 

	my $tryps = return15merPeptide($line[$seqCol], $line[$codeCol], $line[$posCol]);
	if($tryps eq '0'){
	    print "Found error at line #: $lineNum with $line[$codeCol]\t $line[$posCol]\n";
	    }
	$line .= "\t$tryps\n";
	print OUT $line;
	$lineNum +=1;
    
    }
    close(DATA);
    close(OUT);
    
}

1;
