my $log = "hello.txt";
open(LOG, ">$log") or die "Can't write $log";
print LOG "Hello world";
close LOG;
