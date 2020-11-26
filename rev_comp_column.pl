use strict;
use warnings;

die "\nGiven a comma separated .csv file the script will take the the column given in input and transform it into the reverse complement.\n\nUsage: rev_comp_column.pl <input> <column> <flowcell ID> <output>\n\n" unless @ARGV == 4;

my ($input, $selected_column, $flowcell, $output) = @ARGV;

my $counter = 0;
my @read_rows = ();
my @transformed_rows = ();

open(my $IN, "<$input") or die "error open $input for reading";
open(my $OUT, ">$output") or die "error creating $output";

while (my $line = <$IN>) {
	if ($counter > 1 && $line =~ /,/) {
		my @columns = split /,/, $line;
		$columns[$selected_column] = rev_comp($columns[$selected_column]);
		my $new_line = concatenate_row(\@columns);
		print $OUT $new_line;
	}
	elsif ($counter == 0) {
		print $OUT "[Data]\n";
		$counter++;
	}
	elsif ($counter == 1) {
		print $OUT "FCID,Lane,SampleID,SampleRef,index,index2,SampleName,Control,Recipe,Operator,Project\n";
		$counter++;
	}
}

close $OUT;
close $IN;

sub rev_comp {
	my $str_to_revcomp = shift;
	my $rev_comped = reverse $str_to_revcomp;
	$rev_comped =~ tr/ATGCatgc/TACGtacg/;
	return $rev_comped;
}

sub concatenate_row {
	my $row_ref = shift;
	my @row_as_array = @$row_ref;
	my $concatenated_row = $flowcell;
	for (my $i = 1; $i < scalar(@row_as_array); $i++) {
		$concatenated_row = $concatenated_row . "," . $row_as_array[$i];
	}
	return $concatenated_row;
}
