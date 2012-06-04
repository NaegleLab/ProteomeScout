#! /usr/bin/perl -w
use strict;
use pELMTools;
use fileTools;

if (scalar(@ARGV) < 3){
    print "Usage: ./appendpELMFile.pl pELMfile outputFile numberResidues\n";
    exit;
}

my $inputFile = $ARGV[0];
my $outputFile = $ARGV[1];
my $numberResidues = $ARGV[2];

appendAlignedToFile($inputFile, $outputFile, $numberResidues);


# open(IN, $inputFile) || die "Can't open $inputFile\n";
# if(!-e $outputFile) {system("touch $outputFile");}
# open(OUT, ">$outputFile") || die "Can't open $outputFile\n";
# my $line = <IN>; #chomp header

# #print new header to output file
# my $header = returnHeader($inputFile);
# chomp $header;
# $header = $header."\tpep:aligned";
# open(OUT, "> $outputFile") or die "Can't open output file $outputFile";
# print OUT "$header\n";

# #Get column numbers from header and 
# my $CODE_COL = returnColumnNumber($inputFile, "site:code");
# if($CODE_COL < 0){
#     print "ERROR: Your file $inputFile does not have a site:code header type\n";
#     exit;
# }

# #GET header columns
# my $SEQ_COL = returnColumnNumber($inputFile, "sequence:full");
# my $POS_COL = returnColumnNumber($inputFile, "site:pos");
# if($SEQ_COL < 0 or $POS_COL < 0){
#     print "ERROR: Your file $inputFile does not have a  header type for sequence:full or site:pos\n";
#     exit;
# }
# my $ACC_COL = returnColumnNumber($inputFile, "acc:swiss");


# while(defined($line = <IN>)){
#     chomp $line;
#     my @lineArr = split("\t", $line);
#     my $shortSeq = returnAlignedSequence($numberResidues, $lineArr[$SEQ_COL], $lineArr[$CODE_COL], $lineArr[$POS_COL]);
#     print(OUT $line);
#     print(OUT "\t$shortSeq\n");
# }
# close(IN);
# close(OUT);
