#! /usr/bin/perl -w
use strict;

my %multis;
open (IN, "conversion") or die "Can't read conversion";
while(defined(my $line = <IN>)){
    if ($line =~ /\[/){
	my @line = split(" ", $line);
	$multis{$line[0]} = ();
	my @multi = split("", $line[1]);
	foreach my $index (1..$#multi-1){
	    $multis{$line[0]}->{$multi[$index]} = 1;
	}
    }
}
close IN;
open (IN, "$ARGV[0]") or die "Can't open $ARGV[0]";
my @log = <IN>;
my  @newlog;
close IN;

foreach my $i (0..$#log-1){
    # DO that.
    my $match_any = 0;
    my @i = split(" ", $log[$i]);
    my @ipat = split("", $i[1]);
    my $isize = $i[3];
    foreach my $j ($i+1..$#log){
	my @j = split(" " , $log[$j]);
	my @jpat = split("", $j[1]);
	my $jsize = $j[3];
	my $match = 1;
	for my $index (0..$#ipat) {
	    if (!(((($ipat[$index] eq $jpat[$index]) || ($ipat[$index] eq ".") || (defined $multis{$ipat[$index]}->{$jpat[$index]})) && ($isize == $jsize)))){
		$match = 0;
	    }
	}
	if ($match == 1) {
	    $match_any = 1;
	    print STDERR "$log[$i] eliminated by $log[$j]\n";
	    last;
	}
    }
    if ($match_any == 0){
	print STDERR "$log[$i] not eliminated.\n";
	push @newlog, $log[$i];
    }
}

push @newlog, $log[$#log];

my @newlog2;

foreach my $i (reverse(1..$#newlog)){
    # DO that.
    my $match_any = 0;
    my @i = split(" ", $newlog[$i]);
    my @ipat = split("", $i[1]);
    my $isize = $i[3];
    foreach my $j (reverse(0..$i-1)){
	my @j = split(" " , $newlog[$j]);
	my @jpat = split("", $j[1]);
	my $jsize = $j[3];
	my $match = 1;
	for my $index (0..$#ipat) {
	    if (!(((($ipat[$index] eq $jpat[$index]) || ($ipat[$index] eq ".") || (defined $multis{$ipat[$index]}->{$jpat[$index]})) && ($isize == $jsize)))){
		$match = 0;
	    }
	}
	if ($match == 1) {
	    $match_any = 1;
	    print STDERR "$newlog[$i] eliminated by $newlog[$j]\n";
	    last;
	}
    }
    if ($match_any == 0){
	print STDERR "$newlog[$i] not eliminated.\n";
	push @newlog2, $newlog[$i];
    }
}

# BAJ 3/18/08 This used to say (wrongly):
# push @newlog2, $log[0];
push @newlog2, $newlog[0];

print reverse(@newlog2);

