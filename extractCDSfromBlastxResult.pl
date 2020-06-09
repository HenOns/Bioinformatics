#!/usr/bin/perl
# By Henning
use strict;
use warnings;

die "\nReads in the Trinity assembly and then cuts out the coding region, shifts it to the coding strand and directs it in the positive direction.\n\nUsage: extractCDSfromBlastxResult.pl <nucleotide data> <nucleotide vs proteins blastx> <output>\n\n" unless @ARGV == 3;

my ($input, $blastx, $output) = @ARGV;

my %assembly = (); # Keys corresponding to contig headers point to contig sequence in nucleotides
my %header2cds = ();

my $header = "";
my $concatenated_seq = "";

open(my $IN, "<$input") or die "error open $input for reading";

while (my $line = <$IN>) {
    if ($line =~ /^>(\s*\S+).*/) {
		$header = $1;
		$concatenated_seq = "";
    }
    else {
		chomp($line);
		$concatenated_seq = $concatenated_seq . $line;
		$assembly{$header} = $concatenated_seq;
    }
}

close $IN;

open(my $OUT, ">$output") or die "error creating $output";
open($IN, "<$blastx") or die "error open $blastx for reading";

my $counter = 0;

while (my $line = <$IN>) {
	chomp($line);
	my @blastx_arr = split '\t', $line;
	my $temp_qstart = $blastx_arr[6];
	my $temp_qstop = $blastx_arr[7];
	if ($temp_qstop > $temp_qstart) {
		$temp_qstop = $temp_qstop + 3; # To get the stop codon nucleotides
		$temp_qstart = $temp_qstart - 1;
		my $cds = substr $assembly{$blastx_arr[0]}, $temp_qstart, $temp_qstop - $temp_qstart;
		$cds = uc($cds);
		print $OUT ">$blastx_arr[0]_$counter\n";
		print $OUT $cds . "\n";
		$counter++;
	}
	else {
		my $temp_contig = $assembly{$blastx_arr[0]};
		$temp_contig = reverse $temp_contig;
		$temp_contig =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
		my $temp_length = length($temp_contig);
		$temp_qstop = $temp_qstop - 3; # -3 to get the stop codon nucleotides
		$temp_qstop = $temp_length - $temp_qstop + 1;
		$temp_qstart = $temp_length - $temp_qstart;
		my $cds = substr $temp_contig, $temp_qstart, $temp_qstop - $temp_qstart;
		$cds = uc($cds);
		print $OUT ">$blastx_arr[0]_$counter\n";
		print $OUT $cds . "\n";
		$counter++;
	}
}

close $IN;
close $OUT
