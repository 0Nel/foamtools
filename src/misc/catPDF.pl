#!/usr/bin/perl -w

use strict;


my $filename = $ARGV[0];
my @PDFfiles;
my $args;

unless ( `which pdftk` ) {
	print "error, missing pdftk\n";
	exit ( 0 );
}

foreach my $file ( <*> ) {
	if ( -f $file && $file=~ /.pdf/ ) {
		push ( @PDFfiles, $file );	
	}
}

sort ( @PDFfiles );

$args = join ( " ", @PDFfiles );

print "$args \n";

unless ( my $error = `pdftk $args cat output $filename` ) {
	print "$! \n";
}

if ( -e $filename ) {
	print "successfully concatenated all pdfs! \n";
	if ( `rm $args` ) {
		print "could not remove pdf parts..\n";
	}
}
