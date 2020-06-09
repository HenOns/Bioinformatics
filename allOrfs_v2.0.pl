#!/usr/bin/perl
# By Henning
use strict;
use warnings;
use threads;

die "\nTranslates nucleotide sequences in fasta format into all possible reading frames and prints out the CDS as well.\n\nUsage: allOrfs_v2.0.pl <assembly> <peptide length cut off> <codon table> <threads> <output prefix>\n\n" unless @ARGV == 5;

my ($assembly, $lengthCutoff, $codonTable, $threads, $outputPrefix) = @ARGV;

my $cutOff = 80; # Threshold for when an alternative start should be tried if longest methionine starting peptide is shorter than this

my %codon2aa = ();
codonAssignement();
my %fastaHeader2contigSequence = ();
my @thr = ();
my @aminoAcids = ();
my @cds = ();

my $fastaHeader = "";
my $concatenatedSeq = "";
my $threadID = 0;

# To read in the nucleotides contigs into the hash %fastaHeader2contigSequence where fasta headers point to nucleotide seuquences
open(my $IN, "<$assembly") or die "error open $assembly for reading";

while (my $line = <$IN>) {
    if ($line =~ /^>(\s*\S+).*/) {
		$fastaHeader = $1;
		$concatenatedSeq = "";
		$threadID++;
		if ($threadID == $threads) {
			$threadID = 0;
		}
    }
    else {
		chomp($line);
		$concatenatedSeq = $concatenatedSeq . $line;
		$fastaHeader2contigSequence{$threadID}{$fastaHeader} = $concatenatedSeq;
    }
}

close $IN;

for (my $i = 0; $i < $threads; $i++) {
    $thr[$i] = threads->create({'context' => 'list'}, \&translateORFs, $fastaHeader2contigSequence{$i});
}

for (my $i = 0; $i < $threads; $i++) {
    my ($aa_ref, $cds_ref) = $thr[$i]->join();
    push(@aminoAcids, $aa_ref);
    push(@cds, $cds_ref);
}

my $aa_output = "translated_$outputPrefix.fasta";
my $cds_output = "cds_$outputPrefix.fasta";

open(my $OUTAA, ">$aa_output") or die "error creating $aa_output";
open(my $OUTCDS, ">$cds_output") or die "error creating $cds_output";

for (my $i = 0; $i < scalar(@aminoAcids); $i++) {
	foreach my $h (keys %{$aminoAcids[$i]}) {
		if ($lengthCutoff < length($aminoAcids[$i]{$h})) {
			print $OUTAA ">" . $h . "\n";
			print $OUTAA $aminoAcids[$i]{$h} . "\n";
			print $OUTCDS ">" . $h . "\n";
			print $OUTCDS $cds[$i]{$h} . "\n";
		}
	}
}

close $OUTAA;
close $OUTCDS;

sub translateORFs {
	my $fastaHeader2contigSeq_ref = shift;
	my %header2cds = ();
	my %header2aa = ();
	my @readingFrames = (1, 2, 3, -1, -2, -3);
	for (my $i = 0; $i < scalar(@readingFrames); $i++) {
		foreach my $originalHeader (keys %{$fastaHeader2contigSeq_ref}) {
			my $temp_nt = $$fastaHeader2contigSeq_ref{$originalHeader};
			$temp_nt = uc($temp_nt);
			$temp_nt = adjustToFrame($temp_nt, $readingFrames[$i]);
			my $temp_cds = "";
			my $temp_aa = "";
			my $contigID = 1;
			while ($temp_nt =~ /^(.{3})(.*)/) {
				my $temp_codon = $1;
				$temp_nt = $2;
				$temp_cds = $temp_cds . $temp_codon;
				unless (defined $codon2aa{$temp_codon}) {
					$codon2aa{$temp_codon} = "X";
				}
				$temp_aa = $temp_aa . $codon2aa{$temp_codon};
				if ($codon2aa{$temp_codon} eq "*") {
					my ($temp_peptide, $temp_header) = assignPeptide($temp_aa, $originalHeader, $readingFrames[$i], $contigID);
					$contigID++;
					$header2cds{$temp_header} = $temp_cds;
					$header2aa{$temp_header} = $temp_peptide;
					$temp_cds = "";
					$temp_aa = "";
				}
			}
			# To handle the amino acids translated after the last stop or in the case where there are no stops
			$contigID++;
			my $modifiedHeader = "";
			my $peptide = "";
			if ($temp_aa =~ /^[GAVLIFMPWSTYCNQDEHKRXliv]*(m[GAVLIFMPWSTYCNQDEHKRXmliv]*)/) {
				my $end = $1;
				if (length($end) > $cutOff) {
					$peptide = uc($end);
					$modifiedHeader = $originalHeader . "_[Frame_$readingFrames[$i]]_[Stop_unassigned]_$contigID";
				}
				elsif ($temp_aa =~ /^[GAVLIFMPWSTYCNQDEHKRX]*([mliv][GAVLIFMPWSTYCNQDEHKRXmliv]*)/) {
					$end = $1;
					if (length($end) > $cutOff) {
						$peptide = uc($end);
						$modifiedHeader = $originalHeader . "_[Frame_$readingFrames[$i]]_[Stop_unassigned]_$contigID";
					}
					else {
						$peptide = uc($temp_aa);
						$modifiedHeader = $originalHeader . "_[Frame_$readingFrames[$i]]_[Both_unassigned]_$contigID";
					}
				}
				$header2cds{$modifiedHeader} = $temp_cds;
				$header2aa{$modifiedHeader} = $peptide;
			}
			else {
				$peptide = uc($temp_aa);
				$modifiedHeader = $originalHeader . "_[Frame_$readingFrames[$i]]_[Both_unassigned]_$contigID";
				$header2cds{$modifiedHeader} = $temp_cds;
				$header2aa{$modifiedHeader} = $peptide;
			}
		}
	}
	return (\%header2aa, \%header2cds);
}

sub assignPeptide { # Remove the star in later versions
	my ($temp_aa, $originalHeader, $frame, $contigID) = @_;
	my $modifiedHeader = "";
	my $peptide = "";
	if ($temp_aa =~ /^[GAVLIFMPWSTYCNQDEHKRXliv]*(m[GAVLIFMPWSTYCNQDEHKRXmliv]*\*)/) {
		my $end = $1;
		if (length($end) > $cutOff) {
			$peptide = uc($end);
			$modifiedHeader = $originalHeader . "_[Frame_$frame]_$contigID";
		}
		elsif ($temp_aa =~ /^[GAVLIFMPWSTYCNQDEHKRX]*([mliv][GAVLIFMPWSTYCNQDEHKRXmliv]*\*)/) {
			$end = $1;
			if (length($end) > $cutOff) {
				$peptide = uc($end);
				$modifiedHeader = $originalHeader . "_[Frame_$frame]_$contigID";
			}
			else {
				$peptide = uc($temp_aa);
				$modifiedHeader = $originalHeader . "_[Frame_$frame]_[Start_unassigned]_$contigID";
			}
		}
	}
	# Consider adding another elseif statement here to deal with frames that lacks "M"
	else {
		$peptide = uc($temp_aa);
		$modifiedHeader = $originalHeader . "_[Frame_$frame]_[Start_unassigned]_$contigID";
	}
	return ($peptide, $modifiedHeader);
}

sub adjustToFrame {
	my ($temp_nt, $frame) = @_;
	if ($frame < 0) {
		$temp_nt =~ tr/ABCDGHMNRSTUVWXY/TVGHCDKNYSAABWXR/;
		$temp_nt = reverse $temp_nt;
	}
	my $nuclToRemove = abs($frame) - 1;
	$temp_nt =~ s/^.{$nuclToRemove}//s;
	return $temp_nt;
}

sub codonAssignement {
    if ($codonTable == 1) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'l',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => '*',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => '*',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'l',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 4) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'l', 'TTG' => 'l',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => '*',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => 'W',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'l',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'i', 'ATC' => 'i', 'ATA' => 'i',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'v',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 5) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'l',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => '*',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => 'W',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'i', 'ATC' => 'i', 'ATA' => 'm',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'S', 'AGG' => 'S',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'v',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 6) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'L',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => 'Q', 'TAG' => 'Q',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => '*',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 10) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'L',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => '*',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => 'C',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 11) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'l',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => '*',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => '*',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'l',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'i', 'ATC' => 'i', 'ATA' => 'i',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'v',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 15) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'L',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => '*', 'TAG' => 'Q',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => '*',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    elsif ($codonTable == 99) {
		%codon2aa = (
		'TTT' => 'F', 'TTC' => 'F',
		'TTA' => 'L', 'TTG' => 'l',
		'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
		'TAT' => 'Y', 'TAC' => 'Y',
		'TAA' => 'Q', 'TAG' => 'Q',
		'TGT' => 'C', 'TGC' => 'C',
		'TGA' => 'W',
		'TGG' => 'W',
		'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'l',
		'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
		'CAT' => 'H', 'CAC' => 'H',
		'CAA' => 'Q', 'CAG' => 'Q',
		'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
		'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
		'ATG' => 'm',
		'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
		'AAT' => 'N', 'AAC' => 'N',
		'AAA' => 'K', 'AAG' => 'K',
		'AGT' => 'S', 'AGC' => 'S',
		'AGA' => 'R', 'AGG' => 'R',
		'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
		'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',
		'GAT' => 'D', 'GAC' => 'D',
		'GAA' => 'E', 'GAG' => 'E',
		'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G'
		);
    }
    else {
		die "\nError: Table $codonTable is not supported in the script.\n\n"
    }
}
