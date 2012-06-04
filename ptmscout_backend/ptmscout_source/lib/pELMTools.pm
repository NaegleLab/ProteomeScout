use strict;
use entrezTools;
use scansiteParser;

# appendPELMKinaseToFile($pELMFile, $dataFile, $outputFile)
# Append phosphoELM kinase annotation to a text file, must have tryps peptide sequence in a column in data field
# Inputs: $pELMFile - pELM file with appended tryps column (pep:tryps)
#         $dataFile - tab separated data file
#         $outputFile - file to write to with pELM kinase annotation
# Kristen Naegle
# January 16, 2008
sub appendPELMKinaseToFile($$$){
    my $pELMFile = shift;
    my $dataFile = shift;
    my $outputFile = shift;
    
    open(FH_IN, $dataFile) || die "Can't open data file $dataFile for appending phosphoELM info\n";
    if(!-e $outputFile){`touch $outputFile`;}
    open(FH_OUT, ">$outputFile") || die "Can't open output file $outputFile for writing \n";
    my $keyCol = returnColumnNumber($pELMFile, "pep:tryps");
    my $dataCol = returnColumnNumber($dataFile, "pep:tryps");
    my $kinCol = returnColumnNumber($pELMFile, "kin:pELM");
# read in the pELM hash based on pep:tryps
    my $pELMHash = returnPELMHash($pELMFile, $keyCol);

# #do the same for the datafile, but check for miscleavages. 
# #my $pepTryps = "IQPAGNTsPR";
    my @totalKinases;

    # From here for each line in datafile, get annotation col
    my $header = <FH_IN>;
    chomp $header;
    $header .= "\tkin:pelm\n";
    print FH_OUT $header;
    while(defined(my $line = <FH_IN>)){
	chomp $line;
	my $trypsPep = returnField($line, $dataCol);
	$trypsPep =~ s/pY/y/g;
	$trypsPep =~ s/pS/s/g;
	$trypsPep =~ s/pT/t/g;
	
	my $kinaseAnn = returnKinaseAnnotation($pELMHash, $kinCol, $trypsPep);
	$line .= "\t$kinaseAnn\n";
	print FH_OUT $line;
    }
    
    close(FH_IN);
    close(FH_OUT);

}

# $kinase = returnKinaseAnnotation($pELMHash, $kinCol, $pepTryps)
# Returns kinase annotation for a trypsinized peptide 
# Inputs: $pELMHash - hash with keys pep:tryps and value is the line of the pELMHash for the tryps peptide
#         $kinCol - column of pELM file of kinase annotation
#         $pepTryps - trypsinized peptide for which to return the kinase annotation for
# Outputs: $kinase - reference to array of kinases for that peptide
# Kristen Naegle
# January 18, 2008
sub returnKinaseAnnotation($$$){
    my $pELMHash = shift;
    my $kinCol = shift;
    my $pepTryps = shift;
    print "Checking $pepTryps\n";
    my @kin;    
    if($pELMHash->{$pepTryps}){
	print "Found IT\n";
	foreach my $item (@{$pELMHash->{$pepTryps}}){
	    my @item = split("\t", $item);
	    my $kinase = $item[$kinCol];
	    if($kinase =~ /;/){
		my $kinases = returnFields($item, $kinCol);
		push @kin, @$kinases;
	    }
	    elsif(($kinase eq " ") or !$kinase){
		#do nothing
	    }
	    else{
		push @kin, $kinase;
	    }
	}
	if(!@kin){
	   push @kin, "~~~";
	}
	else{
	    my $kin = uniq(\@kin);
	    @kin = @$kin;
	}
	
    }
    else{

	push @kin, "---";
    }
    my $kinField = catArray(\@kin);
    return $kinField;


    

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

# $trypsPeptide = returnTrypsPeptide($sequence, $code, $position);
# Returns the ideal trypsinized peptide around the position and code given in the sequence
# Inputs: $sequence the full AA sequence 
#         $code can be S, T, or Y represents the phosphorylated residue
#         $alignNum is the position in the full sequence (ones-based) of the $code 
# Test file: appendTrypsFragment.pl
# Kristen Naegle, 1/14/08
sub returnTrypsPeptide($$$){
#    my $numAA = 7;
    my $strSequence = shift;
    my $code = shift;
    my $alignNum = shift;

    my @seq = split(//, $strSequence);
    my $pos = $alignNum - 1; #shift to 0-based counting
    
    #if the position and the code are aligned then return the substring
    if($seq[$pos] eq $code){
	my $numAAL = 1;
	while(!($seq[$pos-$numAAL] =~ /[KR]/) && ($pos-$numAAL >= 0)) {
	    $numAAL += 1;
	    
	}
	$numAAL -= 1; #remove the k/r on left
	my $numAAR = 1;
	while(!($seq[$pos+$numAAR] =~ /[KR]/) && ($pos+$numAAR <= $#seq)){
	    $numAAR += 1;
	}
		
	my $posStart = $pos - $numAAL;
	my $posEnd = $pos + $numAAR;
	$seq[$pos] = lc($seq[$pos]); #if phosphorylated change to a lowercase 
	my @subSeq = @seq[$posStart..$posEnd];
	my $subStrSeq = join("", @subSeq);
	
	return $subStrSeq;
    }
    # else if they don't align return zero
    else {
	print "In sub returnAlignedSequences $seq[$pos] does not equal $code\n";
	return 0;
    }
}


# $trypsPeptide = return15merPeptide($sequence, $code, $position);
# Returns the 15mer peptide (+/-7aa around site of phosphorylation)
# Inputs: $sequence the full AA sequence 
#         $code can be S, T, or Y represents the phosphorylated residue
#         $alignNum is the position in the full sequence (ones-based) of the $code 
# Outputs: $trypsPeptide - 15mer peptide around site of phosphorylation (tryps is throwback to old code)
# Kristen Naegle
# Oct. 30, 2008
sub return15merPeptide($$$){
    my $numAA = 7;
    my $strSequence = shift;
    my $code = shift;
    my $alignNum = shift;

    my @seq = split(//, $strSequence);
    my $pos = $alignNum - 1; #shift to 0-based counting
    
    #if the position and the code are aligned then return the substring
    if($seq[$pos] eq $code){
	my $posStart = $pos - $numAA;
	if($posStart <0){
	    $posStart = 0;
	}
	my $posEnd = $pos + $numAA;
	if($posEnd > $#seq){
	    $posEnd = $#seq;
	} 
	$seq[$pos] = lc($seq[$pos]); #if phosphorylated change to a lowercase 
	my @subSeq = @seq[$posStart..$posEnd];
	my $subStrSeq = join("", @subSeq);
	
	return $subStrSeq;
    }
    # else if they don't align return zero
    else {
	print "ERROR: In sub returnAlignedSequences $seq[$pos] does not equal $code\n";
	return 0;
    }
}


# # printKinaseGroup($inputfile, $outputFile, $kinaseCODE, $numAA, $code)
# # prints a SORTED list of aligned sequences (based on $code) of length (2*$numAA + 1) that match kinases 
# # inputs:  $inputFile - tab dileneated file that has the following sequence (no meta labels)
# #                      Gene Symbol	Acc	Sequence	Pos	Code	PubMed ID	Kinase	Expt	Date
# #          $outputFile - file to be written to in OVERWRITE mode, will only write sequences to match up to motif search
# #          $kinaseCODE - string of kinases seperated by commas (example: "EGFR,ErbB2") (Comparison on lowercase forms)
# #          $numAA - length n-terminal and c-terminal side to extend from 
# #          $code - S, Y or T, ensure there's no error in the kinases by doubly checking this 
# # outputs: \@strArray - reference to array of substrings of phosphopeptides targeted by kinase
# # TestFile createEGFRSubgroups.pl 
# # Kristen Naegle, 4/12/07
# # sub printKinaseGroup($$$$$){
# #     my $inputFile = shift;
# #     my $outputFile = shift;
# #     my $kinaseCode = shift;
# #     my $numAA = shift;
# #     my $code = shift;
# #     my $numberSubstrates = 0;
# #     my @strArr; 

# #     #remove whitespaces from any of the kinases and lower case them
# #     $kinaseCode =~ s/ //g;
# #     $kinaseCode = lc($kinaseCode);
# #     my @kinases = split(/,/, $kinaseCode);
    
# #     open(IN, $inputFile) || die "Can't open input file $inputFile\n";
# #     if(!-e $outputFile) { system("touch $outputFile");}
# #     open(OUT, ">$outputFile") || die "Can't open output file for writing $outputFile\n";
# #     my $line;
# #     while(defined($line = <IN>)){
# # 	$line =~ s/ //g; #strip spaces to make sure comparisons work
# # 	$line =~ s/\n//g;
# # 	my @lineArr = split("\t", $line);
# # 	my $lineKinase = $lineArr[KINASE_COL];
# # 	$lineKinase = lc($lineKinase);
# # 	#Don't bother checking if the kinase field is blank
# # 	if($lineKinase){
# # 	    if($code eq $lineArr[CODE_COL]){
# # 		#Wprint "$lineKinase\n";
# # 		foreach my $item (@kinases){
# # 		    if($lineKinase =~ /$item/) {
# # 			my $subString = returnAlignedSequence($numAA, $lineArr[SEQ_COL], $code, $lineArr[POS_COL]);
# # 			push(@strArr, $subString);
# # 			#print "$subString matched in $lineKinase\n";
# # 			$numberSubstrates += 1;
# # 			foreach my $toPrint (@lineArr) {
# # 			    print(OUT "$toPrint\t");
# # 			} 
# # 			print(OUT "\n");
# # 			last;
# # 		    } # end if 
# # 		} #end foreach 
# # 	    } #end check if codes match
# # 	} #end if lineKinase
# # 	chomp $line;
# #     } #end while
    
    
# # #    print("You asked to look for these kinases: @kinases\n");

# #     close(IN);
# #     close(OUT);
# #     @strArr = sort(@strArr);
# #     return \@strArr;

# # }





# # # printAlignedSequencesToFile($inputFile, $outputFile, $numAA, $code)
# # # Takes a tab dilneated input file (appended form) and prints aligned sequences to file for use as background or foreground 
# # # inputs:  $inputFile - tab dileneated file that has the following sequence (no meta labels)
# # #                      Gene Symbol	Acc	Sequence	Pos	Code	PubMed ID	Kinase	Expt	Date
# # #          $outputFile - file to be written to in OVERWRITE mode, will only write sequences to match up to motif search
# # #          $numAA - length n-terminal and c-terminal side to extend from 
# # #          $code - S, Y or T, ensure there's no error in the kinases by doubly checking this 
# # # Test file: createEGFRSubgroups.pl
# # # Kristen Naegle 4/12/07
# # sub printAlignedSequencesToFile($$$$){
# #     my $inputFile = shift;
# #     my $outputFile = shift;
# #     my $numAA = shift;
# #     my $code = shift;
# #     open(IN, $inputFile) || die "Can't open input file $inputFile\n";
# #     my $line;
# #     my @string;
# #     #stuff them into an array 
# #     while(defined($line = <IN>)){
# # 	$line =~ s/ //g; #strip spaces to make sure comparisons work
# # 	$line =~ s/\n//g;
# # 	my @lineArr = split("\t", $line);
# # 	my $subString = returnAlignedSequence($numAA, $lineArr[SEQ_COL], $code, $lineArr[POS_COL]);
# # 	push(@string, $subString);
# #     }
# # #print array to the file
# #     @string = sort(@string);
# #     printArrayToFile($outputFile, \@string);

# # }


# printArrayToFile($outputFile, \@array);
# writes an array to outputFile with new lines after each element
# Inputs: $outputFile - file to write array to
#         $array - reference to array 
# Kristen Naegle
sub printArrayToFile($$){
    my $outputFile = shift;
    my $rArray = shift;
    my @array = @$rArray;
    if(!-e $outputFile) {
	system("touch $outputFile"); 
#	print("Creating new file to print array to $outputFile\n");
    }
    #else { print "Appending array to $outputFile\n";}
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    
    foreach my $item (@array){
	print(OUT "$item\n");
    }
    close(OUT);

}

# $pELMHash = returnPELMHash($pELMFile, $keyCol)
# Returns a hash with key according to $keyCol and value is line of file
# Inputs: $pELMFile - tab separated pELM file
#         $keyCol - column number to use for key (i.e. pep:tryps)
# Outputs: $pELMHash - hash reference 
# Kristen Naegle
# January 18, 2008
sub returnPELMHash($$){
   my $pELMFile = shift;
   my $keyCol = shift;
   my %hash;
   open(FX_67, $pELMFile) || die "Can't open phosphoELM file $pELMFile for reading into hash\n";
   my $line = <FX_67>;
   while(defined($line = <FX_67>)){
       chomp $line;
       my @line = split("\t", $line);
       my $key = $line[$keyCol];
       $key =~ s/ //g;
       if(not $hash{$key}){
	   $hash{$key} = [];
       }
       push @{$hash{$key}}, $line;
       #$hash{$key} = $line;
   }
   

   close(FX_67);
   return \%hash;
}

# $hash = returnPepTrypsHash($file, $keyCol)
# Returns a hash with key according to $keyCol and value is line of file
# Inputs: $file - tab separated pELM file
#         $keyCol - column number to use for key (i.e. pep:tryps)
# Outputs: $hash - hash reference, keys have replaced pY with y, etc.
# Kristen Naegle
# January 18, 2008
sub returnPepTrypsHash($$){
    my $file = shift;
    my $keyCol = shift;

    my %hash;

    open(FX_98, $file) || die "Can't open $file file for reading \n";
    my $line = <FX_98>;
    while(defined($line = <FX_98>)){
	chomp $line;
	my @line = split("\t", $line);
	my $key = returnTrypsPepKey($line[$keyCol]);
	$key =~ s/pY/y/g;
	$key =~ s/pS/s/g;
	$key =~ s/pT/s/g;
	my @keys;
	my ($pos, $codes, $pepCap) = returnPhosphoPos($key);
	if(scalar(@$pos) > 1){
	    foreach my $i (0..scalar(@$pos)-1){
		my @single = split("",$pepCap);		
		$single[$pos->[$i]] = lc($single[$pos->[$i]]);
		my $single = join("", @single);
		#print "Position: ".$pos->[$i]." is one and has form: $single\n";
		push @keys, $single;
	    }
	}
	else{
	    push @keys, $key;
	}
	foreach my $k (@keys){
	    if(not $hash{$k}){
		$hash{$k} = [];
	    }
	    push @{$hash{$k}}, $line;
	}

    }
    return \%hash;
    close(FX_98);

}


# $trypsPepKey = returnTrypsPepKey($pep)
# Returns the idealized trypsinization of a peptide
# Inputs: $pep - peptide to find trypsinization of (phosphorylations denoted by pY, pS, pT or y, s, t)
# Outputs: - $trypsPep - idealized trypsinized peptide
# Kristen Naegle
# January 10, 2008
sub returnTrypsPepKey($){
    my $pep = shift;
    
    $pep =~ s/pY/y/;
    $pep =~ s/pS/s/;
    $pep =~ s/pT/s/;


    #check to see if multiple phosphorylations

    my $foundSubstr = 0;
    my $pepTryps;
    my @sub = split(/[RK]/, $pep);
    if(scalar(@sub) > 1){
	foreach my $i (0..$#sub){
	    my $substr = $sub[$i];
	    if($substr =~ /[a-z]/){
		$foundSubstr = 1;
		my $index = index($pep, $substr);
		my @s = split("", $pep);
		$pepTryps = $substr.$s[$index+length($substr)]; 
		last;
	    }
	}
    }
    else{
	$pepTryps = $pep;
    }
    $pepTryps =~ s/ //g;
    return $pepTryps;

} 


# appendAlignedToFile($inputFile, $outputFile, $numberResidues);
# Opens $inputFile and writes to $outputFile with appended column of the aligned sequence based on pELM input format with sequence:full site:pos site:code
# Inputs: $inputFile - tab separated input file with sequence and positions
#         $outputFile -destination of file
#         $numberResidues - number of residues to align on each side
# Kristen Naegle
# Nov. 30, 2009
sub appendAlignedToFile($$$){
    my ($inputFile, $outputFile, $numberResidues) = @_;
    open(IN, $inputFile) || die "Can't open $inputFile\n";
    if(!-e $outputFile) {system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open $outputFile\n";
    my $line = <IN>; #chomp header
    
#print new header to output file
    my $header = returnHeader($inputFile);
    chomp $header;
    $header = $header."\tpep:aligned";
    open(OUT, "> $outputFile") or die "Can't open output file $outputFile";
    print OUT "$header\n";
    
#Get column numbers from header and 
    my $CODE_COL = returnColumnNumber($inputFile, "site:code");
    if($CODE_COL < 0){
	print "ERROR: Your file $inputFile does not have a site:code header type\n";
	exit;
    }
    
#GET header columns
    my $SEQ_COL = returnColumnNumber($inputFile, "sequence:full");
    my $POS_COL = returnColumnNumber($inputFile, "site:pos");
    if($SEQ_COL < 0 or $POS_COL < 0){
	print "ERROR: Your file $inputFile does not have a  header type for sequence:full or site:pos\n";
	exit;
    }
    my $ACC_COL = returnColumnNumber($inputFile, "acc:swiss");
    
    
    while(defined($line = <IN>)){
	chomp $line;
	my @lineArr = split("\t", $line);
	my $shortSeq = returnAlignedSequence($numberResidues, $lineArr[$SEQ_COL], $lineArr[$CODE_COL], $lineArr[$POS_COL]);
	print(OUT $line);
	print(OUT "\t$shortSeq\n");
    }
    close(IN);
    close(OUT);
    
}

1;
