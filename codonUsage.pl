#!/usr/bin/perl
# codonUsage.pl
# By Henning Onsbring (henning.onsbring@icm.uu.se)
use strict;
use warnings;

die "\n    Script analyze codon usage by identifiyng conserved proteins
    by BLASTX. These proteins are aligned in multiple MAFFT alignments.
    The multiple sequence alignments are then used to score the
    confidence of the start and stop codon usage inferred from BLASTX.\n
    Requirements:
    - Installation of mafft (http://mafft.cbrc.jp/alignment/software/linux.html)
    - Installation of BLAST (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/)
    - Both MAFFT and BLAST in your UNIX path
    - Database in fasta format with conserved proteins
    - Headers in the database formated according to '>IDUniqueForProtein_IDUniqueForOrganism\n
    Usage: codonUsage.pl <Transcriptome assembly> <Conserved proteins database> <Threads> <Outputfolder>\n\n" unless @ARGV == 4;

my ($input, $database, $nthreads, $outputfolder) = @ARGV;

my %codon2aa = ();
my %start_codon2no_ocurrances = ();
my %stop_codon2no_ocurrances = ();
my %non_start_or_stop_codon2no_ocurrances = ();
my %primary_start = ();
my %secondary_starts = ();
my %primary_stop = ();
my %secondary_stops = ();
initiate_codon_statistics();

my %header2bitscore = (); # Fasta header points to the bitscore from the blastx alignment against the database with conserved proteins
my %header2ntseq = (); # Fasta header points to the nucleotide sequence from the trinity assembly
my %protein_id2if_seen = (); # Protein id as key points to if that particular id has been seen by the loop that goes through the sorted blastx result
my %header2if_seen = (); # Fasta header points to if that particular contig has been seen by the loop that goes through the sorted blastx result
my %dbheader2aaseq = (); # Fasta header in the protein database points to the corresponding amino acid sequence
my %header2complement_info = (); # Fasta header points to if the nucleotide sequence should be changed to the complement or not
my %header2start_coord = (); # Fasta header points to start coordinate of the nucleotide sequence that align to the protein database
my %header2stop_coord = (); # Fasta header points to stop coordinate of the nucleotide sequence that align to the protein database
my %header2protein_id = (); # Fasta header points to the protein identifier that it got the best blast hit to
my %header2identity = (); # Fasta header points to the percent identity of the best alignment hit
my %non_redundant_headers2ntseq = (); # Fasta header points to the corresponding sequence. Hash is non redundant for the nucleotide sequences of the protein hits
my %header2coding_ntseq = (); # Fasta header points to the corresponding coding nucleotide sequence. The hash is non reduntant, i.e. only one nucleotide sequence per protein from the database
my %header2start_codon = (); # Fasta header points to what start codon the corresponding sequence use
my %header2stop_codon = (); # Fasta header points to what stop codon the corresponding sequence use
my %header2aa_seq = (); # Fasta header points to the corresponding translated predicted coding nucleotide sequence. The hash is non reduntant, i.e. only one sequence per protein from the database
my %protein_id2path_to_seq_file_for_aln = (); # Protein ID points to the path of a file containing all protein sequences of the corresponding ID that will be used in a multiple sequence alignment
my %header2start_ratio = (); # Fasta header points to the percentage this particular sequence share the start with other sequences in the alignment
my %header2stop_ratio = (); # Fasta header points to the percentage this particular sequence share the stop with other sequences in the alignment

# To add an / at the end out the outputfolder string if that is not given by the user
if ($outputfolder =~ /^.+[^\/]$/) {
    $outputfolder = $outputfolder . "/";
}

# Creates a blast database of the conserved proteins
my $blastdb = $outputfolder . "codon_usage_conserved_proteins.db";
system("makeblastdb -in $database -dbtype prot -out $blastdb");

# Run blastx with the assembeled data as query and the conserved proteins as database
# Table 6 is a good starting guess since the table only use one stop codon that is otherwise tryptophan, which is a rare amino acid
my $pluss_blastx_output = $outputfolder . "assembly_pluss_vs_conserved_proteins_t6.blastx";
my $minuss_blastx_output = $outputfolder . "assembly_minuss_vs_conserved_proteins_t6.blastx";
system("blastx -query $input -db $blastdb -num_threads $nthreads -out $pluss_blastx_output -query_gencode 6 -strand plus -outfmt '6 std salltitles'");
system("blastx -query $input -db $blastdb -num_threads $nthreads -out $minuss_blastx_output -query_gencode 6 -strand minus -outfmt '6 std salltitles'");

# Mark columns with minusstrand and plusstrand and get rid of blasthits with too low significanse by creating a refined blastoutput, > 60% sequence ID and alignment length longer than 70 amino acids
my $p_refined_blastx_output = $outputfolder . "assembly_refined_pluss_vs_conserved_proteins_t6.blastx";
my $m_refined_blastx_output = $outputfolder . "assembly_refined_minuss_vs_conserved_proteins_t6.blastx";
system("cat $pluss_blastx_output | grep -E \"^\\S+\\s+\\S+\\s+100|^\\S+\\s+\\S+\\s+[6-9][0-9]\" | grep -E \"^\\S+\\s+\\S+\\s+\\S+\\s+[7-9][0-9]|^\\S+\\s+\\S+\\s+\\S+\\s+[0-9]{3}\" | sed -E 's/^(\\S+\\s+[^_]+)_\\S+(\\s+.+)/\\1\\2/' | sed -E 's/^(\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+)\\S+/\\1p/' > $p_refined_blastx_output");
system("cat $minuss_blastx_output | grep -E \"^\\S+\\s+\\S+\\s+100|^\\S+\\s+\\S+\\s+[6-9][0-9]\" | grep -E \"^\\S+\\s+\\S+\\s+\\S+\\s+[7-9][0-9]|^\\S+\\s+\\S+\\s+\\S+\\s+[0-9]{3}\" | sed -E 's/^(\\S+\\s+[^_]+)_\\S+(\\s+.+)/\\1\\2/' | sed -E 's/^(\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+)\\S+/\\1m/' > $m_refined_blastx_output");

# Make a list of all the contigs in the assembly that had a significant blasthit to the database
my $concatenated_refined_blastx = $outputfolder . "concatenated_pm_refined.blastx";
system("cat $p_refined_blastx_output $m_refined_blastx_output | sort -n -r -k 12 -k 4 -k 3 -k 2 > $concatenated_refined_blastx");
my $contig_list = $outputfolder . "conserved_contigs.list";
system("cat $concatenated_refined_blastx | cut -f 1 | sort | uniq > $contig_list");

# Create a refined assembly with only the contigs that aligned in the blast search
my $refined_assembly = $outputfolder . "refined_assembly.fasta";

# To create the file refined_assembly.fasta containing only the nucleotide sequences with significant blast hits
my %header2refined_assembly_ntseq = ();
my %header2save_info = ();
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
		$header2refined_assembly_ntseq{$header} = $concatenated_seq;
		$header2save_info{$header} = 0;
    }
}

close $IN;

open($IN, "<$contig_list") or die "error open $contig_list for reading";

while (my $line = <$IN>) {
    chomp($line);
    $header2save_info{$line} = 1;
}

close $IN;

open(my $OUT, ">$refined_assembly") or die "error creating $refined_assembly";

foreach my $h (keys %header2refined_assembly_ntseq) {
    if ($header2save_info{$h} == 1) {
	print $OUT ">" . $h . "\n";
	my $temp_seq = uc($header2refined_assembly_ntseq{$h}) . "\n";
	while ($temp_seq =~ /^(.{60})(.*)/) {
	    print $OUT $1 . "\n";
	    $temp_seq = $2;
	}
		if ($temp_seq =~ /^(.+)$/) {
			print $OUT $1 . "\n";
		}
    }
}

close $OUT;

# To build a hash where the keys, i.e. headers in the fasta file, points at the nucleotide sequence of the assembeled contig
open($IN, "<$refined_assembly") or die "error open $refined_assembly for reading";

while (my $line = <$IN>) {
    if ($line =~ /^>(\s*\S+).*/) {
		$header = $1;
		$concatenated_seq = "";
		$header2bitscore{$header} = 0;
    }
    else {
		chomp($line);
		$concatenated_seq = $concatenated_seq . $line;
		$header2ntseq{$header} = $concatenated_seq;
    }
}

close $IN;

# To build a hash where the keys, i.e.headers in the fasta file, points at the corresponding amino acid sequence of that contig in the database with conserved proteins
open($IN, "<$database") or die "error open $database for reading";

while (my $line = <$IN>) {
    if ($line =~ /^(>.+)/) {
		$header = $1;
		$concatenated_seq = "";
    }
    else {
		chomp($line);
		$concatenated_seq = $concatenated_seq . $line;
		$dbheader2aaseq{$header} = $concatenated_seq;
    }
}

close $IN;

# To find stop and start coordinates of the blasthits on the nucleotide contigs in the concatenated refined assembly. Also this is neede to identify what contigs that needs to be changed to the complement sequence before translated to amin acids.
open($IN, "<$concatenated_refined_blastx") or die "error open $concatenated_refined_blastx for reading";

while (my $line = <$IN>) {
    if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)/) {
		my $temp_header = $1;
		my $temp_hit_protein_id = $2;
		my $temp_identity = $3;
		my $temp_start = $4;
		my $temp_stop = $5;
		my $temp_bitscore = $6;
		my $temp_strand_info = $7;
		if ($temp_start > $temp_stop) {
			$temp_stop = $temp_stop - 3; # To get the stop codon nucleotides
		}
		else {
			$temp_stop = $temp_stop + 3; # To get the stop codon nucleotides
		}
		unless (exists $protein_id2if_seen{$temp_hit_protein_id}) {
			unless (exists $header2if_seen{$temp_header}) { # Need modification if implementing support for genomic data
				if ($temp_strand_info eq "p") {
					$header2complement_info{$temp_header} = "n"; # no need to change to the complement
				}
				elsif ($temp_strand_info eq "m") {
					$header2complement_info{$temp_header} = "c"; # no need to change to the complement
				}
				if ($temp_start > $temp_stop) {
					$header2ntseq{$temp_header} = reverse $header2ntseq{$temp_header};
					my $temp_length = length($header2ntseq{$temp_header});
					$header2start_coord{$temp_header} = $temp_length - $temp_start;
					$header2stop_coord{$temp_header} = $temp_length - $temp_stop + 1;
				}
				else {
					$header2start_coord{$temp_header} = $temp_start - 1;
					$header2stop_coord{$temp_header} = $temp_stop;
				}
				$header2bitscore{$temp_header} = $temp_bitscore;
				$header2protein_id{$temp_header} = $temp_hit_protein_id;
				$header2identity{$temp_header} = $temp_identity;
				$protein_id2if_seen{$temp_hit_protein_id} = "seen"; # If the key exist in the hash the protein id has been seen by the loop before
				$header2if_seen{$temp_header} = "seen"; # If the key exist in the hash the contig has been seen by the loop before
			}
		}
    }
}

close $IN;

# Loop will change contigs identified in the %header2complement_info hash into the complement if needed
foreach my $h (keys %header2ntseq) {
    if (exists $header2complement_info{$h} && $header2complement_info{$h} eq "c") {
		$header2ntseq{$h} =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    }
}

# If multiple contigs in the assembly align to the same protein these lines make sure only the best hit continues to be analyzed
my %protein_id2score_4comparison = ();
my %protein_id2score_best = ();
my %protein_id2header = ();
my %protein_id2seq = ();

foreach my $h (keys %header2protein_id) {
    $protein_id2score_best{$header2protein_id{$h}} = 0; # To allow if the current score is better than the best score observed so far
}

############# Could possibly be removed #############

# $header2protein_id{$h} can several times lead to the same protein id. Iterate over all contigs so the hit to each protein id with best bit score gets stored and overwrite the hits with lower score
foreach my $h (keys %header2protein_id) {
    $protein_id2score_4comparison{$header2protein_id{$h}} = $header2bitscore{$h};
    if ($protein_id2score_4comparison{$header2protein_id{$h}} > $protein_id2score_best{$header2protein_id{$h}}) {
		$protein_id2score_best{$header2protein_id{$h}} = $protein_id2score_4comparison{$header2protein_id{$h}};
		$protein_id2header{$header2protein_id{$h}} = $h;
		$protein_id2seq{$header2protein_id{$h}} = $header2ntseq{$h};
    }
}

foreach my $pid (keys %protein_id2seq) {
    my $temp_seq = $protein_id2seq{$pid};
    my $temp_header = $protein_id2header{$pid};
    $non_redundant_headers2ntseq{$temp_header} = $temp_seq;
}

####################################################

# Loop cuts out the coding sequence and the stop codon
foreach my $h (keys %non_redundant_headers2ntseq) {
    $header2coding_ntseq{$h} = substr $header2ntseq{$h}, $header2start_coord{$h}, $header2stop_coord{$h} - $header2start_coord{$h};
}

# To create an output file with the predicted coding parts of the assembeled nucleotide sequences
my $cds_output = $outputfolder . "cds_output.fasta"; # Path for the file containing the coding sequences from the refined assembly
open($OUT, ">$cds_output") or die "error creating $cds_output";

foreach my $h (keys %header2coding_ntseq) {
    if (length($header2coding_ntseq{$h}) > 27) {
		print $OUT ">" . $h . "\n";
		print $OUT $header2coding_ntseq{$h} . "\n";
    }
    else {
		delete $header2bitscore{$h};
		delete $header2ntseq{$h};
		delete $dbheader2aaseq{$h};
		delete $header2complement_info{$h};
		delete $header2start_coord{$h};
		delete $header2stop_coord{$h};
		delete $header2protein_id{$h};
		delete $header2identity{$h};
		delete $non_redundant_headers2ntseq{$h};
		delete $header2coding_ntseq{$h};
    }
}

close $OUT;

# To count codon frequenzy
foreach my $h (keys %header2coding_ntseq) {
    if ($header2coding_ntseq{$h} =~ /^(.{3})(.+)(.{3})$/) {
		my $temp_start = $1;
		my $temp_coding = $2;
		my $temp_stop = $3;
		$header2start_codon{$h} = $1;
		$header2stop_codon{$h} = $3;
		$start_codon2no_ocurrances{$temp_start}++;
		$stop_codon2no_ocurrances{$temp_stop}++;
		while ($temp_coding =~ /^(.{3})(.*)/) {
			my $temp_codon = $1;
			$temp_coding = $2;
			$non_start_or_stop_codon2no_ocurrances{$temp_codon}++;
			#Can be used to check stop codons that are potentially miss predicted
			#if ($temp_codon eq "TAA" || $temp_codon eq "TAG") {
			#print $head . "\n";
			#}
		}
    }
}

my @stop_codons_sorted = ();
my @start_codons_sorted = ();
my @non_start_or_stop_codons_sorted = ();

foreach my $codon (sort {$stop_codon2no_ocurrances{$b} <=> $stop_codon2no_ocurrances{$a}} keys %stop_codon2no_ocurrances) {
    push (@stop_codons_sorted, $codon . " " . $stop_codon2no_ocurrances{$codon});
}

foreach my $codon (sort {$start_codon2no_ocurrances{$b} <=> $start_codon2no_ocurrances{$a}} keys %start_codon2no_ocurrances) {
    push (@start_codons_sorted, $codon . " " . $start_codon2no_ocurrances{$codon});
}

foreach my $codon (sort {$non_start_or_stop_codon2no_ocurrances{$b} <=> $non_start_or_stop_codon2no_ocurrances{$a}} keys %non_start_or_stop_codon2no_ocurrances) {
    push (@non_start_or_stop_codons_sorted, $codon . " " . $non_start_or_stop_codon2no_ocurrances{$codon});
}

# To find out what codons that is among top 25 in the preliminary prediction of stops and at the same time top 5 least used codon 
my %pred_stop_codon2top25info = ();
my %pred_sense_codon2last5info = ();

for (my $i = 0; $i < 25; $i++) {
    if ($stop_codons_sorted[$i] =~ /^(\w{3})\s/) {
		my $trunc_line = $1;
		$pred_stop_codon2top25info{$trunc_line} = "top25";
    }
}
for (my $i = 0; $i < 5; $i++) {
    if ($non_start_or_stop_codons_sorted[63-$i] =~ /^(\w{3})\s/) {
		my $trunc_line = $1;
		$pred_sense_codon2last5info{$trunc_line} = "last5";
		print "Last5 sense: $trunc_line\n";
    }
}

foreach my $top_stop (keys %pred_stop_codon2top25info) {
    foreach my $last_sense (keys %pred_sense_codon2last5info) {
		if (exists $pred_stop_codon2top25info{$top_stop} && exists $pred_sense_codon2last5info{$last_sense}) {
			if ($top_stop eq $last_sense) {
				$primary_stop{$top_stop} = "s" # This codon stored in this hash is potentially a stop codon
			}
		}
    }
}

my $codon_frequenzy_output = $outputfolder . "preliminary_codon_frequenzy_estimation.txt";
open($OUT, ">$codon_frequenzy_output") or die "error creating $codon_frequenzy_output";

print "\n### Prelimenary codon usage estimation ###\n\n";
print "Stop:\tStart:\tOther:\n";
for (my $i = 0; $i < scalar(@stop_codons_sorted); $i++) {
    print "$stop_codons_sorted[$i]\t$start_codons_sorted[$i]\t$non_start_or_stop_codons_sorted[$i]\n";
    print $OUT "$stop_codons_sorted[$i]\t$start_codons_sorted[$i]\t$non_start_or_stop_codons_sorted[$i]\n";
}
print "\nStatistics also written to $codon_frequenzy_output\n";

close $OUT;
	
foreach my $h (keys %header2coding_ntseq) {
    my $dna_seq = $header2coding_ntseq{$h};
    my $temp_aas = "";
    while ($dna_seq =~ /^(\w{3})(.+)/) {
		$temp_aas = $temp_aas . $codon2aa{$1};
		$dna_seq = $2;
    }
    $header2aa_seq{$h} = $temp_aas;
}

my $aa_output = $outputfolder . "aa_output.fasta";

open($OUT, ">$aa_output") or die "error creating $aa_output";

foreach my $h (keys %header2aa_seq) {
    print $OUT ">" . $h . "\n";
    print $OUT $header2aa_seq{$h} . "\n";
}

close $OUT;

# Files with amino acid seuqences are created, each file will be used to construct a multiple sequence alignment
foreach my $h (keys %header2aa_seq) {

    my $aln_output = $outputfolder . "aln_" . $header2protein_id{$h}; # Path to the file that will be created and can be alined
    open($OUT, ">$aln_output") or die "error creating $aln_output";

    print $OUT ">" . $h . "\n";
    print $OUT $header2aa_seq{$h} . "\n";

	foreach my $ref_header (keys %dbheader2aaseq) {
	    if ($ref_header =~ /^>$header2protein_id{$h}_/) {
			print $OUT $ref_header . "\n";
			print $OUT $dbheader2aaseq{$ref_header} . "\n";
			$protein_id2path_to_seq_file_for_aln{$header2protein_id{$h}} = $aln_output;
	    }
	}

close $OUT;

}

foreach my $protein_id (keys %protein_id2path_to_seq_file_for_aln) {
    my $temp_output = $protein_id2path_to_seq_file_for_aln{$protein_id} . ".fasta";
    system("mafft --thread $nthreads --maxiterate 1000 --localpair $protein_id2path_to_seq_file_for_aln{$protein_id} > $temp_output");
}

foreach my $protein_id (keys %protein_id2path_to_seq_file_for_aln) {

    my $temp_aln = $protein_id2path_to_seq_file_for_aln{$protein_id} . ".fasta";
    my %header2aa_wgaps = (); # Keys corresponding to fasta headers points to the aligned amino acid sequences with gaps %temp_alignement
    my %ngaps_start2no_occurances = (); # Keys corresponding to a certain number of gaps points to the times this number of gaps occur at the start %temp_start_gap_count
    my %ngaps_stop2no_occurances = (); # Keys corresponding to a certain number of gaps points to the times this number of gaps occur at the stop %temp_stop_gap_count
    
    open($IN, "<$temp_aln") or die "error open $temp_aln for reading";

    while (my $line = <$IN>) {
		if ($line =~ /^(>.+)/) {
			$header = $1;
			$concatenated_seq = "";
		}
		else {
			chomp($line);
			$concatenated_seq = $concatenated_seq . $line;
			$header2aa_wgaps{$header} = $concatenated_seq;
		}
    }

    close $IN;
    
    my $total_number_of_contigs = 0;
    my $start_count_of_interest;
    my $stop_count_of_interest;
    my $header_of_interest = $protein_id2header{$protein_id};
    
    foreach my $h (keys %header2aa_wgaps) {
		$total_number_of_contigs++;
		my $temp_start_gaps = "";
		my $temp_stop_gaps = "";
		my $temp_length_start_gaps = 0;
		my $temp_length_stop_gaps = 0;
		if ($header2aa_wgaps{$h} =~ /[^-]*(-+)$/) {
			$temp_stop_gaps = $1;
			$temp_length_stop_gaps = length($temp_stop_gaps);
		}
		if ($header2aa_wgaps{$h} =~ /^(-+)[^-]*/) {
			$temp_start_gaps = $1;
			$temp_length_start_gaps = length($temp_start_gaps);
		}
		if (exists $ngaps_start2no_occurances{$temp_length_start_gaps}) {
			$ngaps_start2no_occurances{$temp_length_start_gaps}++;
		}
		else {
			$ngaps_start2no_occurances{$temp_length_start_gaps} = 1;
		}
		if (exists $ngaps_stop2no_occurances{$temp_length_stop_gaps}) {
			$ngaps_stop2no_occurances{$temp_length_stop_gaps}++;
		}
		else {
			$ngaps_stop2no_occurances{$temp_length_stop_gaps} = 1;
		}
		if ($h =~ /^>$header_of_interest/) { # The exists part should maybe be removed in later versions when possible, could cause errors
			$start_count_of_interest = $temp_length_start_gaps;
			$stop_count_of_interest = $temp_length_stop_gaps;
		}
    }
    my $temp_start_ratio = $ngaps_start2no_occurances{$start_count_of_interest} / $total_number_of_contigs;
    my $temp_stop_ratio = $ngaps_stop2no_occurances{$stop_count_of_interest} / $total_number_of_contigs;
    $header2start_ratio{$header_of_interest} = $temp_start_ratio;
    $header2stop_ratio{$header_of_interest} = $temp_stop_ratio;
    
}

# Initiate a hash that will be a hash of arrays where the fasta header points to an array where each index of the arrays corresponds to different positions along the nucleotide sequence and if that position contain a TAA, TAG or TGA codon
my %header2array_of_startcodoninfo = ();
my %header2array_of_stopcodoninfo = ();
foreach my $h (keys %header2aa_seq) {
    for (my $i = -10; $i < 11; $i++) {
		$header2array_of_startcodoninfo{$h}[$i+10] = "-";
		$header2array_of_stopcodoninfo{$h}[$i+10] = "-";
    }
}

# This loop will check codons +/- 10 positions from the blastx based prediction of the start and stop
for (my $i = -10; $i < 11; $i++) {
    if ($i == 0) {
		goto ENDOFLOOP; # When $i equals to 0 this corresponds to the codons already predicted based on the blastx search
    }
    foreach my $h (keys %header2aa_seq) {
		my $temp_start_coord = $header2start_coord{$h} + (3*$i);
		my $temp_stop_coord = $header2stop_coord{$h} + (3*$i);
		my $temp_length = length($header2ntseq{$h});
		my $temp_start_cut_out_sequence = "";
		my $temp_stop_cut_out_sequence = "";
		my $temp_start_codon = "";
		my $temp_stop_codon = "";
		if ($temp_start_coord > 2) {
			$temp_start_cut_out_sequence = substr $header2ntseq{$h}, $temp_start_coord, $header2stop_coord{$h} - $temp_start_coord;
		}
		if ($temp_stop_coord <= $temp_length) {
			$temp_stop_cut_out_sequence = substr $header2ntseq{$h}, $header2start_coord{$h}, $temp_stop_coord - $header2start_coord{$h};
		}
		if ($temp_start_cut_out_sequence =~ /^(\w{3}).+/) {
			$temp_start_codon = $1;
			$header2array_of_startcodoninfo{$h}[$i+10] = $temp_start_codon;
		}
		if ($temp_stop_cut_out_sequence =~ /.+(\w{3})$/) {
			$temp_stop_codon = $1;
			$header2array_of_stopcodoninfo{$h}[$i+10] = $temp_stop_codon;
		}
    }
    ENDOFLOOP:
}

my $main_result_start = $outputfolder . "main_result_start_codon_analysis";

open($OUT, ">$main_result_start") or die "error creating $main_result_start";

foreach my $h (sort {$header2start_ratio{$b} <=> $header2start_ratio{$a}}  keys %header2aa_seq) {
    my $temp_sum_of_result = ($header2identity{$h} / 100) + $header2start_ratio{$h};
    print $OUT $h . "\t" . $header2protein_id{$h} . "\t" . $header2start_codon{$h} . "\t" . $header2identity{$h} . "\t" . $header2start_ratio{$h} . "\t" . $temp_sum_of_result . "\t" . $header2array_of_startcodoninfo{$h}[0] . "\t" . $header2array_of_startcodoninfo{$h}[1] . "\t" . $header2array_of_startcodoninfo{$h}[2] . "\t" . $header2array_of_startcodoninfo{$h}[3] . "\t" . $header2array_of_startcodoninfo{$h}[4] . "\t" . $header2array_of_startcodoninfo{$h}[5] . "\t" . $header2array_of_startcodoninfo{$h}[6] . "\t" . $header2array_of_startcodoninfo{$h}[7] . "\t" . $header2array_of_startcodoninfo{$h}[8] . "\t" . $header2array_of_startcodoninfo{$h}[9] . "\t" . $header2array_of_startcodoninfo{$h}[11] . "\t" . $header2array_of_startcodoninfo{$h}[12] . "\t" . $header2array_of_startcodoninfo{$h}[13] . "\t" . $header2array_of_startcodoninfo{$h}[14] . "\t" . $header2array_of_startcodoninfo{$h}[15] . "\t" . $header2array_of_startcodoninfo{$h}[16] . "\t" . $header2array_of_startcodoninfo{$h}[17] . "\t" . $header2array_of_startcodoninfo{$h}[18] . "\t" . $header2array_of_startcodoninfo{$h}[19] . "\t" . $header2array_of_startcodoninfo{$h}[20] . "\n";
}

close $OUT;

my $main_result_stop = $outputfolder . "main_result_stop_codon_analysis";

open($OUT, ">$main_result_stop") or die "error creating $main_result_stop";

foreach my $h (sort {$header2stop_ratio{$b} <=> $header2stop_ratio{$a}}  keys %header2aa_seq) {
    my $temp_sum_of_result = ($header2identity{$h} / 100) + $header2stop_ratio{$h};
    print $OUT $h . "\t" . $header2protein_id{$h} . "\t" . $header2stop_codon{$h} . "\t" . $header2identity{$h} . "\t" . $header2stop_ratio{$h} . "\t" . $temp_sum_of_result . "\t" . $header2array_of_stopcodoninfo{$h}[0] . "\t" . $header2array_of_stopcodoninfo{$h}[1] . "\t" . $header2array_of_stopcodoninfo{$h}[2] . "\t" . $header2array_of_stopcodoninfo{$h}[3] . "\t" . $header2array_of_stopcodoninfo{$h}[4] . "\t" . $header2array_of_stopcodoninfo{$h}[5] . "\t" . $header2array_of_stopcodoninfo{$h}[6] . "\t" . $header2array_of_stopcodoninfo{$h}[7] . "\t" . $header2array_of_stopcodoninfo{$h}[8] . "\t" . $header2array_of_stopcodoninfo{$h}[9] . "\t" . $header2array_of_stopcodoninfo{$h}[11] . "\t" . $header2array_of_stopcodoninfo{$h}[12] . "\t" . $header2array_of_stopcodoninfo{$h}[13] . "\t" . $header2array_of_stopcodoninfo{$h}[14] . "\t" . $header2array_of_stopcodoninfo{$h}[15] . "\t" . $header2array_of_stopcodoninfo{$h}[16] . "\t" . $header2array_of_stopcodoninfo{$h}[17] . "\t" . $header2array_of_stopcodoninfo{$h}[18] . "\t" . $header2array_of_stopcodoninfo{$h}[19] . "\t" . $header2array_of_stopcodoninfo{$h}[20] . "\n";
}

close $OUT;

print "\n### Computing refined codon usage estimation ###\n\n";

correctionOfPreliminaryCodonEstimation($main_result_stop,\%header2array_of_stopcodoninfo, \%primary_stop, \%secondary_stops,  "main_result_stop_codon_analysis.refined");
correctionOfPreliminaryCodonEstimation($main_result_start,\%header2array_of_startcodoninfo, \%primary_start, \%secondary_starts,  "main_result_start_codon_analysis.refined");

print "Refined statistics of start codon usage:\n";
system("cat $outputfolder/main_result_start_codon_analysis.refined | sed -E 's/^\\S+\\s+\\S+\\s+(\\S+)/\\1/' | sort | uniq -c | sort -n -r");
print "\nRefined statistics of stop codon usage:\n";
system("cat $outputfolder/main_result_stop_codon_analysis.refined | sed -E 's/^\\S+\\s+\\S+\\s+(\\S+)/\\1/' | sort | uniq -c | sort -n -r");

print "\nTwo lists created listing every protein identified in the blastx alignment and the corresponding codon prediction. The lists are sorted so the predictions done with highest confidence are starting from the top.\n";
print "Stop codon predictions: $main_result_stop\n";
print "Start codon predictions $main_result_start\n\n";
print "Comparisson of the abundance for start, stop and other codons: $codon_frequenzy_output\n\n";

sub correctionOfPreliminaryCodonEstimation {
    my ($main_result, $updownstream_codons_ref, $primary_codons_ref, $secondary_codons_ref, $filename) = @_;
    my %primary_codons = %$primary_codons_ref;
    my %secondary_codons = %$secondary_codons_ref;
    my %updownstream_codons = %$updownstream_codons_ref;
    my @checking_order = (9, 11, 8, 12, 7, 13, 6, 14, 5, 15, 4, 16, 3, 17, 2, 18, 1, 19, 0, 20); # Index 9 corresponds to one step before the predicted start/stop, index 11 corresponds to one step after the predicted start/stop
    my $refined_result = $outputfolder . $filename;
    open($IN, "<$main_result") or die "error open $main_result for reading";
    open($OUT, ">$refined_result") or die "error creating $refined_result";
    
    while (my $line = <$IN>) {
    
		if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)/) {
			my $temp_header = $1;
			my $temp_protein = $2;
			my $initial_codon = $3;
			my $primary_refined_codon = "";
			my $secondary_refined_codon = "";
			my $final_codon;
			
			my $first_counter = 0;
			for (my $i = 0; $i < scalar(@checking_order); $i++) {
			$first_counter++;
			my $temp_codon = $updownstream_codons{$temp_header}[$checking_order[$i]];
			if (exists $primary_codons{$temp_codon}) {
				$primary_refined_codon = $temp_codon;
				goto OUTOFFIRSTLOOP;
			}
			}
			OUTOFFIRSTLOOP:
			my $second_counter = 0;
			for (my $i = 0; $i < scalar(@checking_order); $i++) {
			$second_counter++;
			my $temp_codon = $updownstream_codons{$temp_header}[$checking_order[$i]];
			if (exists $secondary_codons{$temp_codon}) {
				$secondary_refined_codon = $temp_codon;
				goto OUTOFSECONDLOOP;
			}
			}
			OUTOFSECONDLOOP:
			if (exists $primary_codons{$initial_codon}) {
			$final_codon = $initial_codon;
			}
			elsif (exists $primary_codons{$primary_refined_codon} && $first_counter < 9) {
			$final_codon = $primary_refined_codon;
			}
			elsif (exists $secondary_codons{$secondary_refined_codon} && $second_counter < 7) {
			$final_codon = $secondary_refined_codon;
			}
			elsif (exists $primary_codons{$primary_refined_codon} && $first_counter < 15) {
			$final_codon = $primary_refined_codon;
			}
			elsif (exists $secondary_codons{$secondary_refined_codon} && $second_counter < 13) {
			$final_codon = $secondary_refined_codon;
			}
			elsif (exists $primary_codons{$primary_refined_codon}) {
			$final_codon = $primary_refined_codon;
			}
			elsif (exists $secondary_codons{$secondary_refined_codon}) {
			$final_codon = $secondary_refined_codon;
			}
			else {
			$final_codon = $initial_codon;
			}
			
			print $OUT "$temp_header $temp_protein $final_codon\n";
		}
    }
    
    close $OUT;
    close $IN;
}

sub initiate_codon_statistics {

    %start_codon2no_ocurrances = (
    'TTT' => 0, 'TTC' => 0,
    'TTA' => 0, 'TTG' => 0,
    'TCT' => 0, 'TCC' => 0, 'TCA' => 0, 'TCG' => 0,
    'TAT' => 0, 'TAC' => 0,
    'TAA' => 0, 'TAG' => 0,
    'TGT' => 0, 'TGC' => 0,
    'TGA' => 0,
    'TGG' => 0,
    'CTT' => 0, 'CTC' => 0, 'CTA' => 0, 'CTG' => 0,
    'CCT' => 0, 'CCC' => 0, 'CCA' => 0, 'CCG' => 0,
    'CAT' => 0, 'CAC' => 0,
    'CAA' => 0, 'CAG' => 0,
    'CGT' => 0, 'CGC' => 0, 'CGA' => 0, 'CGG' => 0,
    'ATT' => 0, 'ATC' => 0, 'ATA' => 0,
    'ATG' => 0,
    'ACT' => 0, 'ACC' => 0, 'ACA' => 0, 'ACG' => 0,
    'AAT' => 0, 'AAC' => 0,
    'AAA' => 0, 'AAG' => 0,
    'AGT' => 0, 'AGC' => 0,
    'AGA' => 0, 'AGG' => 0,
    'GTT' => 0, 'GTC' => 0, 'GTA' => 0, 'GTG' => 0,
    'GCT' => 0, 'GCC' => 0, 'GCA' => 0, 'GCG' => 0,
    'GAT' => 0, 'GAC' => 0,
    'GAA' => 0, 'GAG' => 0,
    'GGT' => 0, 'GGC' => 0, 'GGA' => 0, 'GGG' => 0
    );
	    
    %stop_codon2no_ocurrances = (
    'TTT' => 0, 'TTC' => 0,
    'TTA' => 0, 'TTG' => 0,
    'TCT' => 0, 'TCC' => 0, 'TCA' => 0, 'TCG' => 0,
    'TAT' => 0, 'TAC' => 0,
    'TAA' => 0, 'TAG' => 0,
    'TGT' => 0, 'TGC' => 0,
    'TGA' => 0,
    'TGG' => 0,
    'CTT' => 0, 'CTC' => 0, 'CTA' => 0, 'CTG' => 0,
    'CCT' => 0, 'CCC' => 0, 'CCA' => 0, 'CCG' => 0,
    'CAT' => 0, 'CAC' => 0,
    'CAA' => 0, 'CAG' => 0,
    'CGT' => 0, 'CGC' => 0, 'CGA' => 0, 'CGG' => 0,
    'ATT' => 0, 'ATC' => 0, 'ATA' => 0,
    'ATG' => 0,
    'ACT' => 0, 'ACC' => 0, 'ACA' => 0, 'ACG' => 0,
    'AAT' => 0, 'AAC' => 0,
    'AAA' => 0, 'AAG' => 0,
    'AGT' => 0, 'AGC' => 0,
    'AGA' => 0, 'AGG' => 0,
    'GTT' => 0, 'GTC' => 0, 'GTA' => 0, 'GTG' => 0,
    'GCT' => 0, 'GCC' => 0, 'GCA' => 0, 'GCG' => 0,
    'GAT' => 0, 'GAC' => 0,
    'GAA' => 0, 'GAG' => 0,
    'GGT' => 0, 'GGC' => 0, 'GGA' => 0, 'GGG' => 0
    );

    %non_start_or_stop_codon2no_ocurrances = (
    'TTT' => 0, 'TTC' => 0,
    'TTA' => 0, 'TTG' => 0,
    'TCT' => 0, 'TCC' => 0, 'TCA' => 0, 'TCG' => 0,
    'TAT' => 0, 'TAC' => 0,
    'TAA' => 0, 'TAG' => 0,
    'TGT' => 0, 'TGC' => 0,
    'TGA' => 0,
    'TGG' => 0,
    'CTT' => 0, 'CTC' => 0, 'CTA' => 0, 'CTG' => 0,
    'CCT' => 0, 'CCC' => 0, 'CCA' => 0, 'CCG' => 0,
    'CAT' => 0, 'CAC' => 0,
    'CAA' => 0, 'CAG' => 0,
    'CGT' => 0, 'CGC' => 0, 'CGA' => 0, 'CGG' => 0,
    'ATT' => 0, 'ATC' => 0, 'ATA' => 0,
    'ATG' => 0,
    'ACT' => 0, 'ACC' => 0, 'ACA' => 0, 'ACG' => 0,
    'AAT' => 0, 'AAC' => 0,
    'AAA' => 0, 'AAG' => 0,
    'AGT' => 0, 'AGC' => 0,
    'AGA' => 0, 'AGG' => 0,
    'GTT' => 0, 'GTC' => 0, 'GTA' => 0, 'GTG' => 0,
    'GCT' => 0, 'GCC' => 0, 'GCA' => 0, 'GCG' => 0,
    'GAT' => 0, 'GAC' => 0,
    'GAA' => 0, 'GAG' => 0,
    'GGT' => 0, 'GGC' => 0, 'GGA' => 0, 'GGG' => 0
    );

    %codon2aa = (
    'TTT' => 'F', 'TTC' => 'F',
    'TTA' => 'L', 'TTG' => 'L',
    'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S',
    'TAT' => 'Y', 'TAC' => 'Y',
    'TAA' => 'Q', 'TAG' => 'Q',
    'TGT' => 'C', 'TGC' => 'C',
    'TGA' => 'W',
    'TGG' => 'W',
    'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
    'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
    'CAT' => 'H', 'CAC' => 'H',
    'CAA' => 'Q', 'CAG' => 'Q',
    'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R',
    'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
    'ATG' => 'M',
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
    
    %primary_start = (
    'ATG' => 'M'
    );
    %secondary_starts = (
    'TTA' => 'M', 'TTG' => 'M', 'CTG' => 'M', 'ATT' => 'M', 'ATC' => 'M', 'ATA' => 'M', 'GTG' => 'M' # Start codons seen in the litterature
    );
    %secondary_stops = (
    'TAA' => 'stop', 'TAG' => 'stop', 'TGA' => 'stop' # Comon stop codons
    );
    
}
