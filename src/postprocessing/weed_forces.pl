#!/usr/bin/perl -w

use strict;

my $start_time = 0;
my $end_time = 10;

my $pattern = "forces";

## die "please provide a pattern for input\n" unless (scalar @ARGV > 0);

while (defined (my $arg = shift @ARGV)) {
    if ($arg =~ m/\-p/){
        $pattern = shift @ARGV;
    } elsif ($arg =~ m/\-st/) {
        $start_time = shift @ARGV;
    } elsif ($arg =~ m/\-et/) {
		$end_time = shift @ARGV;
	}else {
        print "omitting parameter ", shift @ARGV, "\n";
    }
}

print "input pattern set to: ", $pattern, "\n";
print "start time ist set to: ", $start_time, "\n";
print "end time set to: ", $end_time, "\n";

my $output = "forces_clean.dat";
open (OUT, "> $output") or die "could not open output file handle: $!\n";

my @files = glob("./$pattern*");

foreach ( @files ) {
    # print $_, "\n";
    open (IN, "< $_") or die "could not open file: $!\n";
    while (defined (my $line = <IN>)) {
        chomp $line;
        next if ($line =~ m/\#/);
        $line =~ s/\(/ /g;
        $line =~ s/\)/ /g;
        next if (length($line) == 0);
        print OUT $line, " \n";
    }
}

close (IN) or die "could not close filehandle: $!\n";
close (OUT) or die "could not close fildehandle: $!\n"
