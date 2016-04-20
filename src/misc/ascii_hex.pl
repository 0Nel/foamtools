#!/usr/bin/perl

use strict;

sub ascii_to_hex ($)
{
	## Convert each ASCII character to a two-digit hex number.
	(my $str = shift) =~ s/(.|\n)/sprintf("%02lx", ord $1)/eg;
	return $str;
}

sub hex_to_ascii ($)
{
	## Convert each two-digit hex number back to an ASCII character.
	(my $str = shift) =~ s/([a-fA-F0-9]{2})/chr(hex $1)/eg;
	return $str;
}

sub shift ($ $) {
	my ($str, $num) = shift @_;
	my ($str, $num) = shift @_;
		
}

while (1)
{
	print "Enter a string for conversion: ";
	my $str = <STDIN>;
	chomp $str;
	last if !$str || $str =~ /^(?:quit|end)$/i;
	my $h_str = ascii_to_hex $str;
	print "\n\tHex: $h_str\n\n";
	my $a_str = hex_to_ascii $h_str;
	printf("\tASCII: %s\n\n", $a_str eq $str ? $a_str : "error ($a_str)");
}
