#!/usr/bin/env perl

use strict;
use warnings;

# Main loop: reads line-by-line from STDIN (pipe) or input files
while (my $line = <>) {
    chomp $line;

    # Skip VEP header lines starting with #
    next if $line =~ /^#/;

    # Split the VEP line by tabs
    # $F[0] is the ID, $F[3] is Gene, $F[13] is Extra, etc.
    my @fields = split(/\t/, $line);

    my $chrom  = $fields[0];
    my $start  = $fields[1] + 1;

    my $full_id = $fields[2] // "";  # Prevent undef warning if column 2 is missing
    my $strand = (defined $full_id && $full_id =~ /\+/) ? "+" : "-";

    my $extra_info = $fields[7];
    my $end   = extract_tag($extra_info, 'END');
    my $vep = extract_tag($extra_info, 'CSQ');
    
    # 5. Print the final Tab-Separated line
    print join("\t", $chrom, $start, $end, "INS", $strand, $full_id) . ';' . sort_vep_string($vep) . "\n";
}

# --- Subroutines ---

sub extract_tag {
    my ($text, $tag_name) = @_;
    
    # Attempt to match "TAG_NAME=value" up to the next semicolon or end of string
    if ($text =~ /$tag_name=([^;]+)/) {
        return $1;
    }
    
    # Return dot if not found or empty
    return ".";
}

sub sort_vep_string {
    my ($input_string) = @_;
    
    return "" unless defined $input_string;

    # 1. Split by comma
    my @items = split(/,/, $input_string);

    # 2. Sort numerically
    my @sorted = sort {
        # Extract number from A
        my ($num_a) = extract_number($a);
        
        # Extract number from B
        my ($num_b) = extract_number($b);

        # Compare (Ascending)
        $num_a <=> $num_b
    } @items;

    # 3. Join back together
    return join(',', @sorted);
}

sub extract_number {
    my ($str) = @_;
    
    # Split by pipe
    my @cols = split(/\|/, $str);
    
    # Look for the first field that is purely digits
    # This handles both cases: "|ENSG|Type|668|" and "Gene|ENSG|Type|...|"
    foreach my $col (@cols) {
        if ($col =~ /^\d+$/) {
            return $col; 
        }
    }
    
    # RETURN -1 IF NO NUMBER FOUND (Forces these to the front)
    return -1; 
}