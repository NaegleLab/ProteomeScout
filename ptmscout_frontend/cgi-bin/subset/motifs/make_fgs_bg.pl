#! /usr/bin/perl -w
use strict;

my $cutoff = $ARGV[0];

open(IN, "sites_final.csv") or die "Can't read sites_final.csv";
my $line = <IN>;
chomp $line;
my @line = split('\t', $line);
my %header;
foreach my $i (0..$#line){
    $line[$i] =~ s/^\"//;
    $line[$i] =~ s/\"$//;
    $header{$line[$i]} = $i;
}

open(BG, ">background") or die "Can't write background";
open(FG_M, ">foreground_M-DK") or die "Can't write foreground_M-DK";
open(FG_H, ">foreground_H-DK") or die "Can't write foreground_H-DK";
open(FG_SH, ">foreground_SH-DK") or die "Can't write foreground_SH-DK";
my $bg = 0;  my $m = 0;  my $h = 0;  my $sh = 0;

while(defined(my $line = <IN>)){
    chomp $line;
    my @line = split('\t', $line);
    $line[$header{'pep_aligned'}] =~ s/\"//g;
    print BG "$line[$header{'pep_aligned'}]\n";
    $bg++;
    if ($line[$header{'M/DK'}] >= $cutoff) { $m++; print FG_M "$line[$header{'pep_aligned'}]\n"; }
    if ($line[$header{'H/DK'}] >= $cutoff) { $h++; print FG_H "$line[$header{'pep_aligned'}]\n"; }
    if ($line[$header{'SH/DK'}] >= $cutoff) { $sh++; print FG_SH "$line[$header{'pep_aligned'}]\n"; }
}
close BG; 
close FG_M;
close FG_H;
close FG_SH;
print "BG:\t$bg\nM:\t$m\nH:\t$h\nSH:\t$sh\n";
