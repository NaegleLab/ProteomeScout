use strict;
use warnings;
use fileTools;


# appendScansiteInfo($inputFile, $outputFile, $stringency)
# Writes scansite predictions to file Type:Name:Percentage
# Inputs: $inputFile - tab separated file with pep_aligned field
#         $outputFile - destination of input+scansite predictions
#         $stringency - LOW, MEDIUM, or HIGH, cutoff percentage for scansite
# Kristen Naegle
# July, 2007
sub appendScansiteInfo($$$){
    my $inputFile = shift;
    my $outputFile = shift;
    my $stringency = shift;
    if(!-e $outputFile){system('touch $outputFile');}
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    open(IN, "$inputFile") || die "Can't open $inputFile for reading\n";
#Need head for pep:aligned 
    my $PEP_COL = returnColumnNumber($inputFile, "pep:aligned");


#first append header with scansite header and write to output
    my $header = returnHeader($inputFile);

    chomp $header;
    $header = $header."\tkin:scansite\tbind:scansite\n";
    print(OUT "$header"); #includes new line

    my $line =<IN>; #skip the first line since we just reprinted the new header
    while($line = <IN>){
	chomp $line;
	#my @line = split("\t", $line);
	#my $pep = $line[$PEP_COL];
	my $peps = returnFields($line, $PEP_COL);
	my @peps = @$peps;
	my @kin;
	my @bind;
	my @hashRefArr;
	foreach my $pep (@peps){
	    print "PEPTIDE $pep\n";
	    my $hash = returnScansiteHash($pep, $stringency);
	    push @hashRefArr, $hash;
	}
	
	my $appendage = createScansiteInfoLineMultiPep(\@hashRefArr);
	print "kinase\tbind\n";
	print "$appendage\n";
	print OUT "$line\t$appendage\n";
	
    }
    
    
    close(IN);
    close(OUT);
    
    
}



# getScansiteFromPeptide($pep);
# retrieve and write a scansite hmtl file given an input sequence - writes to fixed file scansiteTemp.htm in current directory
# Inputs: $pep - aligned peptide...requires 15mer
# Outputs: No direct outputs..but writes html code to scansiteTemp.htm use returnScansiteHash to access this info
# Kristen Naegle
# July 11, 2007
sub getScansiteFromPeptide($$){
    my $pep = shift;
    my $stringency = shift;
    my $URL = 'http://scansite.mit.edu/cgi-bin/motifscan_seq?protein_id=TEST&motif_option=all&submit=Submit+Request&domain_flag=0&sequence=*'.$pep.'&stringency='.$stringency;
    my $command = "wget ".'"'.$URL.'"'." -O scansiteTemp.htm";
    print "COMMAND: $command\n";
#    system("wget $URL -O scansiteTemp.htm");
    system("$command");
}

# getScansiteMotifFromPeptide($motifFile, $motifName, $pep, $stringency)
# writes url wget to scansiteTemp.htm for parsing, assumes you've already loaded the motif matrix, and know the location
# Inputs: $motifFile - portion of scansite url that refers to location of motif file
#         $motifName - name of motif you loaded
#         $pep - the 15-mer sequence (oriented at position 8)
#         $stringency - LOW, MEDIUM, or HIGH
# Outputs: Writes scansiteTemp.htm
# Kristen Naegle
# April 4, 2008
sub getScansiteMotifFromPeptide($$$$){
    my ($motifFile, $motifName, $pep, $stringency) = @_;
    my $URL = 'http://scansite.mit.edu/cgi-bin/motifscan_mat?protein_id=13261&motif_option=&submit=Submit+Request&sequence=*'.$pep.'&singlemotif=1&usermotifs=1&mfile1='.$motifFile.'&mname1='.$motifName.'&stringency='.$stringency;
    my $command = "wget ".'"'.$URL.'"'." -O scansiteTemp.htm";
    print "COMMAND: $command\n";
#    system("wget $URL -O scansiteTemp.htm");
    system("$command");
    
}

# $name = parseName($line);
#given an html line return the kinase, sh2 domain protein, or ptb domain protein.  Choose the Gene Card name
# inputs: $line - input line from html file that contains the name of the domain protein or kinase 
# outputs: $name - in gene card format
# Kristen Naegle
# July 11, 2007
sub parseName($){
    my $line = shift;
    my @match; 
    while($line =~ m/<b>(.+?)<\/b>/gi){
	push @match, $1;
    }
    my $name = $match[1]; #Gene Card name 
  #  my @name = split(" ", $name);
    $name =~ s/ //g;
   # $name = $name[0];
    return $name;

}




# ($siteNum, $percentile) = parseSiteLine($line);
# returns the site number that corresponds with scansite alignment and the percentile score for that hit
# inputs: $line - html line that contains site and scores
# outputs: $siteNum - just the number, should be 8 in all cases for 15mer passed to scansite URL
#          $percentile - percentile score. 
# Kristen Naegle
# July 11, 2007
sub parseSiteLine($){
    my $line = shift;
    my $site;
    my @td;
    while ($line =~ m/<td>(.+?)<\/td>/gi){
	push @td, $1;
    }
    $site = $td[0];
    my $percentile = $td[2];
    $percentile =~ s/ //g;
    $percentile =~ s/%//g;
    my $SA = $td[3];
    $site =~ m/(\d)/g;
    my $siteNum = $1;
    $siteNum =~ s/ //g;

    return ($siteNum, $percentile);
    

}



# \%hash = returnScansiteHashForDB($pep);
# Given an aligned peptide (residue 8 is central position), return hash of predicted scansite values
# Inputs: $pep - aligned peptide
# Outputs: \%hash - keys are a count, with array values: (source, value, score);
# hardcoded for a 15mer alignment on the central residue
# Kristen Naegle
# July 11, 2007
sub returnScansiteHashForDB($){
    my ($pep) = @_;
    my $stringency = 'LOW';
    my $scansiteFile = "scansiteTemp.htm";
    my $KINASE = "[K|k]in";
    my $SH2 = "SH2";
    my $PTB = "PTB";
    
    my $SITE = 8;
    
    getScansiteFromPeptide($pep, $stringency);
    
    open(FH_1, $scansiteFile) || die "Can't open file, $scansiteFile, for reading \n";

    my $line;
    my %hash;
    my $source = 'scansite';
    my $count=0;
    
    while (defined($line=<FH_1>)){
	# Check for group line
	#options - bgcolor=red
	while ($line =~ /bgcolor=red/){ #make this a while line 
	    # FOUND Group 
	#    print "FOUND A HEADER\n";
	    my ($group, $type) = getGroup($line);
	#    print "Found group: $group of $type\n";
	    $line = <FH_1>;
	    while(!($line =~ /bgcolor=red/ or $line =~ /table>/)){
		#get geneCard name and percentile for all members of group
		
		#Get GENECARD Name
		
		my $kinaseName = parseName($line);
		#print "Found: $kinaseName\n";
		
		$line = <FH_1>; #SACRIFICIAL Line
		
		#GET Percentage and site
		$line = <FH_1>; #get percent and check 
		my ($siteNum, $percentile) = parseSiteLine($line);
		
		if($siteNum==8){
		    my $sourceT = $source.'_'.$type;
		    my @arr = ($sourceT, $kinaseName."_".$group, $percentile);
		    push @{$hash{$count}}, @arr;
#		    my $t = $kinaseName.":".$group.":".$percentile;
	#	    print "Adding to hash with $type --> $t\n";
		    
		#    push @{$hash{$type}}, $t;
		    $count +=1;
		}
		$line = <FH_1>;
	    } # end parse group

	    
	    # $type gives key for hash
	}
	
	
	# If Group.. process using global FH_1 and stuffing into a hash

    
    }
    close(FH_1);

    if(!$count){
	print "FOUND NO predictons\n";
	my @arr = ($source, '~~~', 'NULL');
	push @{$hash{$count}}, @arr;
    }
    return \%hash;   
}

# \%hash = returnScansiteHash($inputFile);
# Given a scansite .htm file, parse it for the kinases, ptb and sh2 domains. 
# Inputs: $inputFile - scansite html file
# Outputs: \%hash - hashref with up to 3 keys 'kinase', 'sh2' and 'ptb' those point to an array of proteinnames:percentile
# hardcoded for a 15mer alignment on the central residue
# Kristen Naegle
# July 11, 2007
sub returnScansiteHash($$){
    #my $scansiteFile = shift;
    my $pep = shift;
    my $stringency = shift;
    my $scansiteFile = "scansiteTemp.htm";
    my $KINASE = "[K|k]in";
    my $SH2 = "SH2";
    my $PTB = "PTB";
    my $SITE = 8;
    getScansiteFromPeptide($pep, $stringency);
    
    open(FH_1, $scansiteFile) || die "Can't open file, $scansiteFile, for reading \n";

    my $line;
    my %hash;
    #$hash{'kinase'} = [];
    #$hash{'bind'} = [];
    
    while (defined($line=<FH_1>)){
	# Check for group line
	#options - bgcolor=red
	while ($line =~ /bgcolor=red/){ #make this a while line 
	    # FOUND Group 
	#    print "FOUND A HEADER\n";
	    (my $group, my $type) = getGroup($line);
	#    print "Found group: $group of $type\n";
	    $line = <FH_1>;
	    while(!($line =~ /bgcolor=red/ or $line =~ /table>/)){
		#get geneCard name and percentile for all members of group
		
		#Get GENECARD Name
		#$line = <FH_1>; 
		my $kinaseName = parseName($line);
		#print "Found: $kinaseName\n";
		
		$line = <FH_1>; #SACRIFICIAL Line
		
		#GET Percentage and site
		$line = <FH_1>; #get percent and check 
		(my $siteNum, my $percentile) = parseSiteLine($line);
		
		if($siteNum==8){
		    my $t = $kinaseName.":".$group.":".$percentile;
		    print "Adding to hash with $type --> $t\n";
		    
		    push @{$hash{$type}}, $t;
		}
		$line = <FH_1>;
	    } # end parse group

	    
	    # $type gives key for hash
	}
	
	
	# If Group.. process using global FH_1 and stuffing into a hash

	chomp $line;
    
    }
    close(FH_1);
    return \%hash;   
}


# ($group, $type) = getGrou($scansiteGroupLine)
# Parses line to get group and type (kinase or bind)
# Inputs: $line - group line
# Outputs: $group - specific name of predicted group
#          $type - kinase or bind
# Kristen Naegle
# July 2007
sub getGroup($){
    my $line = shift;

    # First parse the group name out of there
    my @b;
    while ($line =~ m/<b>(.+?)<\/b>/gi){
	push @b, $1;
    }
    my @group;
    while ($b[0] =~ m/\((.+?)\)/gi){
	push @group, $1;
    }

    my $group = $group[0];
    my $type;
    if ($group =~ /_[K|k]in/){
	$type = 'kinase';
	
    }
    else{
	$type = 'bind';
    }
	

    return ($group[0], $type);


}

# ### CHANGED 8/8/07 to be more general..this is the original #### 
# # # \%hash = returnScansiteHash($inputFile);
# # # Given a scansite .htm file, parse it for the kinases, ptb and sh2 domains. 
# # # Inputs: $inputFile - scansite html file
# # # Outputs: \%hash - hashref with up to 3 keys 'kinase', 'sh2' and 'ptb' those point to an array of proteinnames:percentile
# # # hardcoded for a 15mer alignment on the central residue
# # # Kristen Naegle
# # # July 11, 2007
# # sub returnScansiteHash($){
# #     #my $scansiteFile = shift;
# #     my $pep = shift;
# #     my $scansiteFile = "scansiteTemp.htm";
# #     my $KINASE = "[K|k]in";
# #     my $SH2 = "SH2";
# #     my $PTB = "PTB";
# #     my $SITE = 8;
# #     getScansiteFromPeptide($pep);
    
# #     open(FH_1, $scansiteFile) || die "Can't open file, $scansiteFile, for reading \n";

# #     my $line;
# #     my %hash;
    
# #     while (defined($line=<FH_1>)){
# # 	if ($line =~ /$KINASE/){
# # 	    print $line."\n\n";
# # 	    chomp $line;
# # 	    $line = <FH_1>;
# # 	    chomp $line; #get to the real first kinase
# # 	    #print "\n\nSHOULD BE KINASE LINE: $line";
# # 	    while($line =~ /$KINASE/){
# # 		my $kinaseName = parseName($line);
# # 		$line = <FH_1>;
# # 		chomp $line;
# # 		#print "\n\nSACRIFICIAL: $line";
# # 		$line = <FH_1>;
# # 		chomp $line;
# # 		#print "\n\n SHOULD BE Y8 LINE: $line\n";
# # 		(my $siteNum, my $percentile) = parseSiteLine($line);
# # 		### ADD ERROR HANDLING 
# # 		if($siteNum eq $SITE){
# # 		    push @{$hash{'kinase'}}, $kinaseName.":".$percentile;
# # 		}
# # 		$line = <FH_1>;
# # 		chomp $line;
# # 	    }
# # 	}
# # 	if ($line =~ /$SH2/){
# # 	    print $line;
# # 	    chomp $line;
# # 	    $line = <FH_1>; 
# # 	    #print "\n\nSHOULD BE SH2 DOMAIN LINE: $line";
# # 	    while($line =~ /$SH2/){
# # 		my $sh2Name = parseName($line);
# # 		chomp $line;
# # 		$line = <FH_1>;
# # 		#print "\n\n SACRIFICIAL LINE: $line";
# # 		$line = <FH_1>;
# # 		#print "\n\n SHOULD BE Y8 Line: $line\n";
# # 		(my $siteNum, my $percentile) = parseSiteLine($line);
# # 		if($siteNum eq $SITE){
# # 		    push @{$hash{'sh2'}}, $sh2Name.":".$percentile;
# # 		}
# # 		$line = <FH_1>;
# # 		chomp $line;
		
		
# # 	    }
	    
# # 	}
# # 	if ($line =~ /$PTB/){
# # 	    print $line;
# # 	    chomp $line;
# # 	    $line = <FH_1>; 
# # 	    #print "\n\nSHOULD BE PTB DOMAIN LINE: $line";
# # 	    while($line =~ /$PTB/){
# # 		my $name = parseName($line);
# # 		chomp $line;
# # 		$line = <FH_1>;
# # 		#print "\n\n SACRIFICIAL LINE: $line";
# # 		$line = <FH_1>;
# # 		#print "\n\n SHOULD BE Y8 Line: $line\n";
# # 		(my $siteNum, my $percentile) = parseSiteLine($line);
# # 		if($siteNum eq $SITE){
# # 		    push @{$hash{'ptb'}}, $name.":".$percentile;
# # 		}
# # 		$line = <FH_1>;
# # 		chomp $line;
		
		
# # 	    }
	    
# # 	}
	
# # 	chomp $line;
    
# #     }
# #     close(FH_1);
# #     return \%hash;   
# #}

# ##CHANGED 8/08/07 to be able to handle parsing scansite more generally
# # # $string = createScansiteInfoLine($scansiteHashRef);
# # # This takes a scansite hash and creates a line that can be appended to a text file
# # # Inputs: $scansiteHashRef - hash reference to scansite object, as created in returnScansiteHash
# # #         keys: 'kinase', 'sh2', 'ptb' **need to create a new function for S/T sites
# # # Outputs: $string - looks like: name:percentile, name2:percentile with tab separation between fields (kinase, sh2 and ptb)
# # # Kristen Naegle
# # # July 12, 2007
# # sub createScansiteInfoLine($){
# #     my $hashRef = shift;
# #     my $kinaseStr;
# #     my $sh2Str;
# #     my $ptbStr;
# #     my $string;
# #     if($hashRef->{'kinase'}){
# # 	$kinaseStr = catArray($hashRef->{'kinase'});
# #     }
# #     else{ $kinaseStr = " "; }
# #     if($hashRef->{'sh2'}){
# #         $sh2Str = catArray($hashRef->{'sh2'});
# #     }
# #     else{ $sh2Str = " ";}
# #     if($hashRef->{'ptb'}){
# # 	$ptbStr = catArray($hashRef->{'ptb'});
# #     }
# #     else{ $ptbStr = " "; }
# #     $string = $kinaseStr."\t".$sh2Str."\t".$ptbStr;
# #     return $string;
    
# # }

# $string = createScansiteInfoLine($scansiteHashRef);
# This takes a scansite hash and creates a line that can be appended to a text file
# Inputs: $scansiteHashRef - hash reference to scansite object, as created in returnScansiteHash
#         keys: 'kinase', 'sh2', 'ptb' **need to create a new function for S/T sites
# Outputs: $string - looks like: name:percentile, name2:percentile with tab separation between fields (kinase, sh2 and ptb)
# Kristen Naegle
# July 12, 2007
sub createScansiteInfoLine($){
    my $hashRef = shift;
    my $kinaseStr;
    my $bindStr; 
    my $string;
    if($hashRef->{'kinase'}){
	$kinaseStr = catArray($hashRef->{'kinase'});
    }
    else{ $kinaseStr = "~~~"; }
    if($hashRef->{'bind'}){
        $bindStr = catArray($hashRef->{'bind'});
    }
    else{ $bindStr = "~~~";}
    $string = $kinaseStr."\t".$bindStr;
    return $string;
}

# Modified create ScansiteInfoLine above, to handle the multiple peptides in a column solution. So call this for each peptide (with corresponding hash) and then create line using blah
# Kristen Naegle
# October 31, 2007
sub returnScansiteInfo($){
    my $hashRef = shift;
    my $kinaseStr;
    my $bindStr; 
    my $string;
    if($hashRef->{'kinase'}){
	$kinaseStr = catArray($hashRef->{'kinase'});
    }
    else{ $kinaseStr = "~~~"; }
    if($hashRef->{'bind'}){
        $bindStr = catArray($hashRef->{'bind'});
    }
    else{ $bindStr = "~~~";}
    #$string = $kinaseStr."\t".$bindStr;
    return ($kinaseStr, $bindStr);

}

# $scansiteLine = createScansiteInforLineMultiPep($hashRefArr)
# Creates scansite info line given multiple peptides
# Inputs: $hashArr - reference to hash array of kinases and binding predictions
# Outputs: $scansiteLine - string for printing to append file
# Kristen Naegle
# June 12, 2007
sub createScansiteInfoLineMultiPep($){
    my $hashRefArr = shift;
    my $kinaseLine = "";
    my $bindLine = "";
    my $count=0;
    foreach my $hashPep (@$hashRefArr){
	$count += 1;
	if ($count > 1){
	    $kinaseLine .= ",";
	    $bindLine .= ",";
	}
	if($hashPep->{'kinase'}){
	   
	    $kinaseLine .= catArray($hashPep->{'kinase'});
	}
	else { $kinaseLine .= "~~~"; } 
	if($hashPep->{'bind'}){
	    $bindLine .= catArray($hashPep->{'bind'});
	}
	else { $bindLine .= "~~~"; } 


    
    }
    my $appendage = $kinaseLine."\t";
    $appendage .= $bindLine;
    return $appendage;

}

# $newStr = removeAboveCutoff($string, $cutoff);
# Takes a string (like what was printed to a file using createScansiteInfoLine) and removes entries (comma seperated) that are above $cutoff
# Inputs: $string - comma seperated line from file (scansite column specifically)
#         $cutoff - percentile cutoff
# Outputs: $newString - $string minus those entries above $cutoff
# Kristen Naegle
# July 12, 2007
sub removeAboveCutoff($$){
    my $str = shift;
    if(!$str){
	return "";
    }
    my $cutoff = shift;
    #split the string first on commas into elements
    my @str = split(";", $str);
    if(not @str){
	push @str, $str;
    }
    my @newStr; 
    #split each element into percentile to determine if it should be thrown away
    foreach my $item (@str){
	(my $name, my $percentile) = split(":", $item);
	if($percentile <= $cutoff){
	    push @newStr, $item;
	}
	
    }
    my $newStr = catArray(\@newStr);
    return $newStr;
}

# $str = catArray($arrayRef);
# create a comma seperated line from an array
# Inputs: $arraRef - reference to array of elements to be concatenated into single line
# Outputs: $str - said string 
# Kristen Naegle
# July 12, 2007
sub catArray($){
    my $arrayRef = shift;
    my $str = $arrayRef->[0];
    my $i;
    for($i=1; $i<scalar(@$arrayRef); $i++){
	$str = $str."; ".$arrayRef->[$i];
	
    }
    return $str;
}

1;
