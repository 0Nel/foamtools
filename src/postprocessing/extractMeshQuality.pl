#!/usr/bin/perl -w

use strict;

use constant TRUE => 1;
use constant FALSE => 0;

### usefull for copy and paste
### (m/^\-(?: |\- )$/)      for parsing parameters

### PROTOTYPES
sub warn($);
sub error($);
sub notify($);
sub help;

### color support for better output
my ($GREEN, $YELLOW, $RED, $NORMAL) = "";
if ( -t STDERR && -t STDOUT) {
    $RED = "[01;31m";
    $YELLOW = "[01;33m";
    $GREEN = "[01;32m";
    $NORMAL = "[0m";
}

my $in_file = "";

if (scalar @ARGV > 0 and -f $ARGV[0]) {
    $in_file = shift @ARGV;
}

my @non_ortho_mean = ();
my @non_ortho_max = ();
my @skew = ();
my @time = ();

open (IN, "< $in_file") or die "could not open file $in_file: $!\n";

while ( defined( my $line = <IN> ) ) {

    if ( $line =~ m/.*Time = (\d+\.\d+).*/) {
        if (defined $1) {
            print "time: ", $1, "\n";
            push (@time, $1);
        }
    }


    if ($line =~ m/.*Mesh non-orthogonality Max: (\d+\.\d+) average: (\d+\.\d+).*/) {
        if (defined $1 and defined $2) {
            print "non ortho max: ", $1, " mean: ", $2, "\n";
            push (@non_ortho_max, $1);
            push (@non_ortho_mean, $2);
        }
    }

    if ($line =~ m/.*Max skewness = (\d+\.\d+).*/) {
        if (defined $1) {
            print "schiefe: ", $1, "\n";
            push (@skew, $1);
        }
    }
}

close (IN) or die "could not close file: $in_file: $! \n";

open (OUT, "> ./out.dat") or die "could not open output file: $!\n";

print OUT "# TIME non-orthogonality_max non-orthogonality_mean skewness\n";

my $counter = 0;

for my $time_step (sort @time) {
    print OUT $time_step, " ", $non_ortho_max[$counter], " ", $non_ortho_mean[$counter], " ", $skew[$counter], "\n";
    $counter ++;
}

close (OUT) or die "could not close file : $!\n";


####################### ENHANCED OUTPUT SUBROUTINES ########################

### used to print formated warnings
sub warn ($) {
    print $RED, "$0: subroutine warn expects an argument!\n" unless (defined ($_[0]));
    print $YELLOW, "$0 $_[0]\n", $NORMAL;
}

### used to print formated error messages: not a substitute to die, but rather an extension
sub error($) {
    print $RED, "$0: subroutine error expects an argument!\n" unless (defined ($_[0]));
    print $YELLOW, "$0 $_[0]\n", $NORMAL;
    exit -1;
}

### used to notify the user; intended to ease debugging
sub notify($) {
    print $RED, "$0: subroutine notify expects an argument!\n" unless (defined ($_[0]));
    print $GREEN, "$0 $_[0]\n", $NORMAL;
}

### the *MANDATORY* help subroutine
sub help {
    print "$0 HELP PAGE\n\nthis is going to be the help screen\n";
    exit 0;
}
