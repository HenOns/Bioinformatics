use strict;
use warnings;

die "\nGiven a comma separated .csv file the script will take the the column given in input and transform it into the reverse complement.\n\nUsage: rev_comp_column.pl <input> <column> <output>\n\n" unless @ARGV == 3;

my ($input, $selected_column, $output) = @ARGV;

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
	else {
		print $OUT $line;
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
	my $concatenated_row = $row_as_array[0];
	for (my $i = 1; $i < scalar(@row_as_array); $i++) {
		$concatenated_row = $concatenated_row . "," . $row_as_array[$i];
	}
	return $concatenated_row;
}
