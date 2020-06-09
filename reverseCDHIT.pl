#!/usr/bin/perl
# By Henning
use strict;
use warnings;

die "\nReads in CD-HIT .clstr file and make a list of all contigs where more than <threshold> contigs are similar per cluster.\n\nUsage: reverseCDHIT.pl <input> <threshold>\n\n" unless @ARGV == 2;

my ($input, $threshold) = @ARGV;

my %cell2cluster = ();
my %cluster2longestcontig = ();

my $cluster = "";
my $length = "";
my $cell = "";
my $contig = "";
my $counter = 0;

open(my $IN, "<$input") or die "error open $input for reading";

while (my $line = <$IN>) {
    if ($line =~ /^(>.+)/) {
	$cluster = $1;
	if ($counter >= $threshold) {
	    foreach my $c (keys %cluster2longestcontig) {
		print $cluster2longestcontig{$c} . "\n";
	    }
	}
	$counter = 0;
	%cell2cluster = ();
	%cluster2longestcontig = ();
    }
    elsif ($line =~ /^\S+\s+\d+nt\S+\s+>([^_]+)_\S+\.\.\.\s+at/) {
	$cell = $1;
	unless (exists $cell2cluster{$cell}) {
	    $counter++;
	}
	$cell2cluster{$cell} = $cluster;
    }
    elsif ($line =~ /^\S+\s+(\d+)nt\S+\s+>([^_]+)(_\S+)\.\.\.\s+\*/) {
	$length = $1;
	$cell = $2;
	my $pre_contig = $3;
	$contig = $cell . $pre_contig;
	unless (exists $cell2cluster{$cell}) {
	    $counter++;
	}
	$cell2cluster{$cell} = $cluster;
	$cluster2longestcontig{$cluster} = $contig;
    }
}

close $IN;
