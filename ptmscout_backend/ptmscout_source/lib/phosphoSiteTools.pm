use strict;
use warnings;

BEGIN{
    use constant    PS_NAME_COL => 0;
    use constant    PS_SWISS_COL => 1;
    use constant    PS_CODE_COL => 2;
    use constant    PS_SPECIES_COL => 3;
    use constant    PS_PEP_COL => 4;
    use constant    PS_HQ_COL => 5;
    use constant    PS_LQ_COL => 6;
    use constant    PS_LIT_COL => 7;
}

# printBkgdnFromPSFile($inputFile
# Prints the sequences from an input file that is columnated like the PhosphoSite data files
# Removes _ and replaces with spaces.  
# Inputs: $inputFile 
# Outputs: writes otuputfile
# Caution..uses phosphosite Hard coded columns
sub printBkgndFromPSFile($){
    my $inputFile = shift;
    open(IN, $inputFile) || die "Can't open $inputFile for reading\n";
    my $outputFile = $inputFile."_bkgnd";
    if(!-e $outputFile){system("touch $outputFile");}
    open(OUT, ">$outputFile") || die "Can't open $outputFile for writing\n";
    my $line;
    while(defined($line = <IN>)){
	my @line = split("\t", $line);
	my $pep = $line[PS_PEP_COL];
	$pep =~ s/_/ /g;
	print(OUT "$pep\n");
	chomp $line;

    }
}

1;
