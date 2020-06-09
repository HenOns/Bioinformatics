#!/usr/bin/perl
# By Henning
use strict;
use warnings;

die "\nGiven a list of fasta headers the script selects these sequences and prints them to a separate file.\n\nUsage: selectSequencesFromFasta.pl <assembly> <list> <output>\n\n" unless @ARGV == 3;

my ($input, $list, $output) = @ARGV;

my %assembly = ();
my $header = "";
my $concatenated_seq = "";
my %save = ();

# To build a hash where the keys (headers in fasta file) points at the assembeled contigs
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
		$save{$header} = 0;
    }
}

close $IN;

open($IN, "<$list") or die "error open $list for reading";

while (my $line = <$IN>) {
    chomp($line);
    $save{$line} = 1;
}

close $IN;

open(my $OUT, ">$output") or die "error creating $output";

foreach my $cheader (keys %assembly) {
    if ($save{$cheader} == 1) {
		print $OUT ">" . $cheader . "\n";
		my $seq = uc($assembly{$cheader}) . "\n";
		while ($seq =~ /^(.{60})(.*)/) {
			print $OUT $1 . "\n";
			$seq = $2;
		}
		if ($seq =~ /^(.+)$/) {
			print $OUT $1 . "\n";
		}
    }
}

close $OUT;
