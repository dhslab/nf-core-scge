#!/usr/bin/env perl
use strict;
use warnings;

# This script reformats VEP tabular output to be compatible with extract_variant_reads_ML.py
# 1. Splits Location into chromosome, start, end
# 2. Moves Uploaded_variation (which contains source info) to the last column
# 3. Prints a compatible header using the provided argument as the last column name

my $last_col_header = $ARGV[0] or die "Usage: $0 <last_column_header>\n";

# Print Header
# extract_variant_reads_ML.py expects: #chromosome, start, end, ... [KeyString]
print "#chromosome\tstart\tend\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra\t$last_col_header\n";

while (<STDIN>) {
    # Skip VEP headers
    next if /^#/;
    chomp;
    
    my @cols = split /\t/;
    
    # Safety check for column count
    if (scalar @cols < 2) {
        next;
    }

    my $info = $cols[0];
    my $loc = $cols[1];
    
    # Parse Location: chr:start-end or chr:start
    my ($chr, $pos_str) = split /:/, $loc;
    my ($start, $end);
    
    if (defined $pos_str && $pos_str =~ /-/) {
        ($start, $end) = split /-/, $pos_str;
    } elsif (defined $pos_str) {
        $start = $pos_str;
        $end = $pos_str;
    } else {
        # Fallback if location format is unexpected
        $chr = $loc;
        $start = 0;
        $end = 0;
    }
    
    # Output: chr, start, end, [All cols except first], Info
    print "$chr\t$start\t$end\t";
    
    # Print cols 1 to end (skipping Uploaded_variation which is at index 0)
    if (scalar @cols > 1) {
        print join("\t", @cols[1..$#cols]);
    }
    
    # Print Info at the end
    print "\t$info\n";
}
