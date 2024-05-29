#!/usr/bin/perl
use strict;
use warnings;

my ($input, $output, $filter) = @ARGV;

my $seq = "";
my $header = "";
open (O, ">$output") or die $!;
if ($input =~ /\.gz$/){
	open (F, "zcat $input |") or die $!;
}
else{
	open (F, "<$input") or die $!;
}
while (<F>){
	if ($_ =~ /^>/){
		print O "$header$seq" if (length($seq) > 0 && $header =~ /$filter/);
		$header = $_;
		$seq = "";
	}
	else{
		$_ =~ tr/U/T/;
		$seq .= $_;
	}
}

print O "$header$seq" if (length($seq) > 0 && $header =~ /$filter/);
exit;
