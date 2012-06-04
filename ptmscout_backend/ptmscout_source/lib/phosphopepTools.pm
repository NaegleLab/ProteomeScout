use strict;
use warnings;
use entrezTools;
use errorHandling;

# ($errorCode, \@phospho) = returnSinglyPhosphoArray($pep)
# from a singly to multiply-phosphorylated peptide, return an array of singly phosphorylated peptides
# Inputs: $pep - where phosphorylations marked by lower case p (preceeding aa) or by lc aa itself
# Outputs: $errorCode - returns true (1) if no phosphorylations found in peptide
#          \@phsopho - array of singly phosphorylated peptides
# Kristen Naegle
# March 9, 2008
sub returnSinglyPhosphoArray($){
    my ($pep) = @_;
    my $errorCode = 0;
    my @phospho;

    $pep =~ s/pT/t/g;
    $pep =~ s/pS/s/g;
    $pep =~ s/pY/y/g;
    if($pep =~ m/[ykst]/g){
	my ($posArr, $code, $pepCap) = returnPhosphoPos($pep);
	for (my $i=0; $i<scalar(@$posArr); $i++){
	    my $p = substr($pepCap, 0, $posArr->[$i]);
	    $p .= lc($code->[$i]);
	    $p .= substr($pepCap, $posArr->[$i]+1, length($pep));
	   # push @phospho, substr($pepCap, $posArr->[$i], length($pep), lc($code->[$i]));
	    push @phospho, $p;
	}
    }
    else{
	$errorCode = 1;
	handleError('returnSinglyPhosphoArray', 'No phosphorylation site', \@_);
    }

    
    return($errorCode, \@phospho);
}

# $trypsPepKey = returnTrypsPep($pep)
# Returns the idealized trypsinization of a peptide
# Inputs: $pep - peptide to find trypsinization of (phosphorylations denoted by pY, pS, pT or y, s, t)
# Outputs: - $trypsPep - idealized trypsinized peptide
# Kristen Naegle
# January 10, 2008
sub returnTrypsPep($){
    my ($pep) = @_;
    chomp $pep;
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




# ($errorCode, $seqPos, $type, $aligned) = returnAlignedandCodeForSinglyPhospho($seq, $pep, $numAA)
# inputs: $seq - peptide sequence of protein
#         $pep - MS peptide fragment (with pS, pY or pT marking alignment)
#         $numAA - number of aa on each side to return
# outputs: $errorCode  1: no phosphorylation in peptide, 2:more than one phosphorylation, 3:peptide doesn't exist in sequence
#          $seqPos - position (aa) in sequence where phosphorylation occurs
#          $type - the amino acid phosphorylated
#          $aligned - corresponding aligned sequences where pS is s pY is y and pT is t
# Kristen Naegle
# March 10, 2008
sub returnAlignedandCodeForSinglyPhospho($$$){
    my ($seq, $pep, $numAA) = @_;

    my $errorCode = 0;

# #now find the position of the phosphorylations 
# #given a peptide in MS format, return capitalized peptide and array of positions
    my $codes;
    my $pos;
    my $pepCap;
    my $aligned;
    my $seqPos;
    my $type;
    ($pos, $codes, $pepCap) = returnPhosphoPos($pep);
    
    if($pos eq "NA"){ #peptide was not found in sequence
	$errorCode = 1; # no phosphorylation site found
	handleError('returnAlignedandCodeForSinglyPhospho', 'Cannot find phosphorylation in peptide', \@_);
    }
    elsif(scalar(@$pos) > 1){
	$errorCode = 2; #more than one phosphorylation site

    }
    else{
	my $index = index($seq, $pepCap);

	if($index == -1){
	    handleError('returnAlignedandCodeForSinglyPhospho', 'Cannot find peptide in sequence', \@_);
	    $errorCode = 1;
	}
	else{
	    my @pos = @$pos;
	    $seqPos = $pos->[0]+$index+1;
	    $aligned = returnAlignedSequence($numAA,$seq, $codes->[0], $seqPos);
	    $type = $codes->[0];
	}
	    
	
    } #end else
    return($errorCode, $seqPos, $type, $aligned);
}



1;
