#!/usr/bin/perl -w

### THIS FILE WAS CREATED BY ALJOSCHA SANDER
### aljoscha.sander@gmail.com
### all content is strictly under GPL v3 


use strict;

sub usage ();
sub get_average (\@ \@ $);
sub average (@);
sub get_stdev (\@ \@ $);
sub stdev (@);
sub print_file (\@ $);

use constant X => 0;
use constant Y => 1;

my @files;
my @data;
my @results;
my $output_file_prefix;
my $n_file = 0;
my $n_line = 0;

if (@ARGV == 0 or $ARGV[0] =~ m/(?:\-h|\-help)/) {
    &usage();
}

while (defined ($_ = shift @ARGV)) {
    if (-f) {
        push (@files, $_);
    } elsif (-d) {
        print "not yet supported\n";
        exit;
    } else {
        $output_file_prefix = $_;
    }
}

foreach my $file (@files) {

    open (FILE, "< $file") or die "could not open file: $!\n";
    while (defined (my $line = <FILE>)) {
        chomp $line;
        my @tmp = split (" ", $line);
        $data[X][$n_line][$n_file] = $tmp[3];
        $data[Y][$n_line][$n_file] = $tmp[4];
        $n_line ++;
    }
    close (FILE) or die "could not close file handle: $! \n";
    $n_line = 0;
    $n_file ++;
}

&get_average(\@data, \@results, $n_file);
&print_file(\@results, "average.txt");

############# SUBROUTINES ###################

sub get_average (\@ \@ $) {
    my ($data, $results, $n_files) = @_;
    my $line = 0;
    while (defined ($data->[X][$line])) {
        $results->[X][$line] = average($data->[X][$line]);
        $results->[Y][$line] = average($data->[Y][$line]);
        $line ++;
    }
}

sub average (@) {
    my @data = @{$_[0]};
    my $length = @data;
    my $sum = 0;
    foreach my $value (@data) {
        $sum += $value;
    }
    return $sum / $length;
}

sub print_file(\@ $) {
    my ($results, $p_file) = @_;
    my $n_line = 0;
    open (OUT, "> $p_file") or die "could not open output file $p_file: $!\n";
    print OUT "X      Y\n";
    while (defined ($results->[X][$n_line]) && defined $results->[Y][$n_line]) {
        print OUT ($results->[X][$n_line]);
        print OUT "    ";
        print OUT ($results->[Y][$n_line]);
        print OUT "\n";
        $n_line++;
    }
    close OUT or die "could not close file: $!\n";
}

sub get_stdev (\@ \@ $) {

}

sub stdev (@) {
    
}

sub usage () {

    print <<EOF;

$0 calculates the average of given files, line by line, position by position
the output file must be specified, otherwise the default file name 'average'
will be used

written by Aljoscha N. Sander, published under the GPLv3
aljoscha.sander\@gmail.com
    
EOF
    exit;
}

