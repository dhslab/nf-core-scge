#!/usr/bin/env python3

from __future__ import division
import edlib
import argparse, re, string, csv, sys
import scipy.stats as stats
import joblib
import pandas as pd
import pyranges as pr
import numpy as np
import pysam
import bisect
from collections import defaultdict

class GenomicDistanceIndex:
    def __init__(self, df):
        """
        Input: Pandas DataFrame with columns 'Chromosome' and 'Pos'
        """
        self.index = {}
        
        # 1. Group by Chromosome and extract positions
        # This is vectorized and significantly faster than iterating rows
        for chrom, group in df.groupby('Chromosome'):
            # Convert the 'Pos' column to a sorted list
            # Sorting is mandatory for binary search (bisect) to work
            self.index[chrom] = sorted(group['Pos'].tolist())

    def get_min_dist(self, query_chrom, query_pos):
        """
        Returns the minimum distance from query_pos to any target on query_chrom.
        Returns float('inf') if chromosome not found.
        """
        # Fast Dictionary Lookup
        if query_chrom not in self.index:
            return float('inf')
        
        positions = self.index[query_chrom]
        
        # Binary Search (O(log N))
        idx = bisect.bisect_left(positions, query_pos)
        
        # Edge Case 1: Insertion point is 0 (closest is the first item)
        if idx == 0:
            return abs(positions[0] - query_pos)
        
        # Edge Case 2: Insertion point is end (closest is the last item)
        if idx == len(positions):
            return abs(positions[-1] - query_pos)
        
        # Normal Case: Check neighbors (left and right of insertion point)
        before = positions[idx - 1]
        after = positions[idx]
        
        return min(abs(before - query_pos), abs(after - query_pos))
    
# =============================================================================
# Constants & Helpers
# =============================================================================

BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_N = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8

CONSUMES_REF = {BAM_CMATCH, BAM_CDEL, BAM_CREF_N, BAM_CEQUAL, BAM_CDIFF}
CONSUMES_READ = {BAM_CMATCH, BAM_CINS, BAM_CSOFT_CLIP, BAM_CEQUAL, BAM_CDIFF}

# Define globals once to avoid rebuilding them on every function call
CIGAR_OPS = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':0, 'X':0}
CIGAR_REGEX = re.compile(r'(\d+)([MIDNSHP=X])')

BND_REGEX = re.compile(r"([ACGTNacgtn]*)(\[|\])([^:]+:\d+)(\[|\])([ACGTNacgtn]*)")

def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

def reverse_complement(seq):
    """Returns the reverse complement of a DNA string."""
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]

def get_sequence(fasta_handle, chrom, start, end):
    """Safe fetch from fasta."""
    if not fasta_handle:
        return None
    try:
        return fasta_handle.fetch(chrom, start, end).upper()
    except (ValueError, KeyError, IndexError):
        return None
    
def parse_sa_tag(sa_string):
    """
    Parses a single SA tag entry.
    Format: rname,pos,strand,CIGAR,mapQ,NM;
    Example: chr2,10500,+,50S100M,60,0;
    """
    parts = sa_string.split(',')
    if len(parts) < 6: return None
    
    return {
        'chrom': parts[0],
        'pos': int(parts[1]), # SA tag positions are 1-based
        'strand': parts[2],   # + or -
        'cigar': parts[3],
        'mapq': int(parts[4]),
        'nm': int(parts[5])
    }

def parse_cigar_string(cigar_str):
    """Parses CIGAR string into list of tuples."""
    if not cigar_str or cigar_str == '*': 
        return []
    # Use findall (which runs in C) and a list comprehension
    return [(CIGAR_OPS[op], int(length)) for length, op in CIGAR_REGEX.findall(cigar_str)]

def cigar_to_string(cigar_tuples):
    codes = "MIDNSHP=X"
    return "".join(f"{l}{codes[op]}" for op, l in cigar_tuples)

def cigar_summary(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_sum = { 'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0 }
    for length, operation in cigar_tuples:
        cigar_sum[operation] += int(length)

    return cigar_sum

def cigar_to_coords(cigar, r_name, r_start_pos, strand, is_primary=False):
    """
    Calculates query (read) and reference coordinates for an alignment.
    Expects CIGAR to be normalized to the Query 5'->3' orientation.
    """
    q_start = 0
    if cigar and cigar[0][0] in [BAM_CSOFT_CLIP, BAM_CHARD_CLIP]:
        q_start = cigar[0][1]
    
    q_consumed_aligned = sum(l for op, l in cigar if op in {BAM_CMATCH, BAM_CINS, BAM_CEQUAL, BAM_CDIFF})
    r_aligned = sum(l for op, l in cigar if op in CONSUMES_REF)
    
    q_end = q_start + q_consumed_aligned - 1

    return {
        'q_start': q_start,
        'q_end': q_end,
        'r_name': r_name,
        'r_start': r_start_pos, 
        'r_end': r_start_pos + r_aligned - 1,
        'strand': strand,
        'is_primary': is_primary
    }

def make_info_string(data):
    """
    Converts a single dict or a list of dicts into a VCF INFO string.
    Aggregates values for shared keys across a list.
    Handles pysam types (tuples for lists, booleans for flags).
    """
    if not data:
        return "."

    # Normalize input: ensure we always iterate over a list of dicts
    if isinstance(data, dict):
        data_list = [data]
    elif isinstance(data, list):
        data_list = data
    else:
        return "."

    info_parts = []
    
    # 1. Identify all unique keys present across the input
    all_keys = sorted(set().union(*(d.keys() for d in data_list)))

    for k in all_keys:
        vals = []
        is_flag = False

        for d in data_list:
            if k in d:
                v = d[k]
                # Handle pysam Flags (boolean True implies presence)
                if isinstance(v, bool):
                    if v: 
                        is_flag = True
                # Handle pysam Tuples/Lists (e.g., AF=0.1,0.2)
                elif isinstance(v, (tuple, list)):
                    vals.extend([str(x) for x in v])
                # Handle scalars (int, str, float)
                else:
                    vals.append(str(v))

        # 2. Format string based on content
        # If it was a flag and we collected no value data (pure flag)
        if is_flag and not vals:
            info_parts.append(k)
        # If we have values, join them (e.g. KEY=val1,val2)
        elif vals:
            info_parts.append(f"{k}={','.join(vals)}")

    return ';'.join(info_parts)

def format_bnd(ref, alt, chrom_mate, pos_mate, fragment, strand_self, strand_mate):
    """
    Formats VCF breakend string.
    strand_self: Strand of the anchor (current VCF record)
    strand_mate: Strand of the target (mate)
    """
    
    # REF ALT Meaning
    # s t[p[ piece extending to the right of p is joined after t
    # s t]p] reverse comp piece extending left of p is joined after t
    # s ]p]t piece extending to the left of p is joined before t
    # s [p[t reverse comp piece extending right of p is joined before t
         
    # Simplified logic matching typical caller output:
    if fragment == "left":
        if strand_self == '+' and strand_mate == '+':
            return f"{ref}{alt}[{chrom_mate}:{pos_mate + 1}["
        elif strand_self == '+' and strand_mate == '-':
            return f"{ref}{alt}]{chrom_mate}:{pos_mate + 1}]"
        elif strand_self == '-' and strand_mate == '-':
            return f"{ref}{alt}[{chrom_mate}:{pos_mate + 1}["
        elif strand_self == '-' and strand_mate == '+':
            return f"[{chrom_mate}:{pos_mate + 1}[{alt}{ref}"
    
    else: # right
        if strand_self == '+' and strand_mate == '+':
            return f"]{chrom_mate}:{pos_mate + 1}]{alt}{ref}"
        elif strand_self == '+' and strand_mate == '-':
            return f"[{chrom_mate}:{pos_mate + 1}[{alt}{ref}"
        elif strand_self == '-' and strand_mate == '-':
            return f"]{chrom_mate}:{pos_mate + 1}]{alt}{ref}"
        elif strand_self == '-' and strand_mate == '+':
            return f"{ref}{alt}]{chrom_mate}:{pos_mate + 1}]"
    
    return f"<{chrom_mate}:{pos_mate + 1}>"

def get_bnd_parts(alt_str):
    """
    Parses VCF BND strings (e.g., ]chr1:123]T) into components.
    Returns: (pre_bases, bracket_char, remote_chrom, remote_pos, post_bases)
    """
    # Regex captures: 1=Pre, 2=Bracket, 3=RemoteLoc, 4=Bracket, 5=Post
    match = BND_REGEX.fullmatch(alt_str)
    if not match:
        return None
    return match.groups()

def generate_contig(chrom,pos,ref_base,alt_base,svtype,ref_file,flank_length):
    """
    Constructs the theoretical alternate sequence (Contig) using the Reference genome.
    Logic:
      1. Identify SV type (Linear vs BND).
      2. Fetch Local Flank (Left or Right depending on break orientation).
      3. Fetch Remote Flank (and RC if necessary based on bracket direction).
      4. Stitch them together.
    """
    # chrom = record['chrom'] #record.chrom
    # pos = record['pos'] #record.pos # 1-based
    # ref_base = record['ref'] #record.ref
    # alt_base = record['alt'] #record.alts[0]
    # svtype = record['type'] #record.info.get('SVTYPE', 'Unknown')
    #strand = record['info']['STRAND'] #record.info.get('STRAND', '')
    
    # Pysam 0-based conversion
    start_idx = pos - 1 
    end_idx = start_idx + len(ref_base)

    try:
        # --- Handle Breakends (BND) ---
        if svtype == 'BND' or '[' in alt_base or ']' in alt_base:
            parts = get_bnd_parts(alt_base)
            if not parts:
                return "Error: Unparseable BND string"

            pre_bases, bracket, remote_loc, _, post_bases = parts
            r_chrom, r_pos = remote_loc.split(':')
            r_pos = int(r_pos)
            
            # 1. Get Local Sequence
            if pre_bases: 
                # Local is upstream (Left) -> Break
                local_seq = ref_file.fetch(chrom, max(0, start_idx - flank_length), end_idx - 1)

                if bracket == '[':
                    r_start_0 = r_pos - 1
                    remote_seq = ref_file.fetch(r_chrom, r_start_0, r_start_0 + flank_length)
                else: # bracket == ']'
                    # Remote points Reverse (Left/Upstream of r_pos) -> RC needed
                    r_end_0 = r_pos
                    raw_remote = ref_file.fetch(r_chrom, max(0, r_end_0 - flank_length), r_end_0)
                    remote_seq = reverse_complement(raw_remote)
                    
                return local_seq + pre_bases + remote_seq #if strand=="+" else reverse_complement(local_seq + pre_bases + remote_seq)

            else: # joined before t (ref is downstream)
                if bracket == '[':
                    local_seq = ref_file.fetch(chrom, start_idx, end_idx + flank_length)
                    # Remote points Forward (Right/Downstream of r_pos)
                    r_start_0 = r_pos - 1
                    raw_remote = ref_file.fetch(r_chrom, r_start_0 + 1, r_start_0 + flank_length)
                    remote_seq = reverse_complement(raw_remote)

                else: # bracket == ']'
                    local_seq = ref_file.fetch(chrom, start_idx + 1, end_idx + flank_length)
                    # Remote points Reverse (Left/Upstream of r_pos) -> RC needed
                    r_end_0 = r_pos
                    remote_seq = ref_file.fetch(r_chrom, max(0, r_end_0 - flank_length), r_end_0)

                return remote_seq + post_bases + local_seq #if strand=="+" else reverse_complement(remote_seq + post_bases + local_seq)

        # --- Handle Linear SVs (DEL, INS, etc) ---
        else:
            
            left_flank = ''
            right_flank = ''
            
            if len(alt_base) > len(ref_base) or svtype == 'INS' or svtype == 'DUP':
                left_flank = ref_file.fetch(chrom, max(0, start_idx - flank_length), start_idx)
                right_flank = ref_file.fetch(chrom, end_idx, end_idx + flank_length)
                # For linear SVs, we sandwich the ALT string between the flanks
            
            elif len(alt_base) < len(ref_base) or svtype == 'DEL':
                left_flank = ref_file.fetch(chrom, max(0, start_idx - flank_length), start_idx)
                right_flank = ref_file.fetch(chrom, end_idx, end_idx + flank_length)
                # For linear SVs, we sandwich the ALT string between the flanks

            seq = f"{left_flank}{alt_base}{right_flank}" # if strand=="+" else reverse_complement(f"{left_flank}{alt_base}{right_flank}")
            
            return seq 

    except KeyError as e:
        return f"Error: Chromosome/Region not found in reference ({e})"
    except Exception as e:
        return f"Error: {str(e)}"

# =============================================================================
# Main function to call SV from a split read alignment. 
# Args are a dict of split alignment coordinates and returns a dict for VCF record generation
# =============================================================================

def call_sv_from_split_read(aln1, aln2, query_seq, fasta_handle=None):
    """
    Core SV calling logic with Ref/Alt sequence generation and Primary alignment anchoring.
    """
    # 1. Sort alignments by their position in the READ (Query Order)
    if aln1['q_start'] < aln2['q_start']:
        L, R = aln1, aln2
    else:
        L, R = aln2, aln1

    read_start = L['q_start']
    read_end = R['q_end']
    read_strand = L['strand'] if L['is_primary'] else R['strand']
    
    # 2. Adjust for Overlap
    read_gap = R['q_start'] - L['q_end'] - 1
    overlap = -read_gap if read_gap < 0 else 0
    
    # Effective Ref coordinates
    # R matches starting at R['r_start'].
    # If there was an overlap, the 'true' breakpoint on the R side is pushed forward/backward.
    if overlap > 0:
        if L['strand'] == R['strand'] or R['strand'] == '+': # if event is a del/dup or a BND and right strand is FWD
            R['r_start'] += overlap
            R['q_start'] += overlap
        else: # if the right strand is - then move it to the left and increase the R start pos.
            R['r_end'] -= overlap
            R['q_start'] += overlap

        read_gap = 0

    # 3. Get SV positions, which are strand dependent. Note ++ and -- are equivalent here.
    bp1_chrom = L['r_name']
    bp1_pos = L['r_end'] if L['strand'] == R['strand'] or L['strand'] == '+' else L['r_start'] # End of first segment (0-based, strand-specific)
    
    bp2_chrom = R['r_name']
    bp2_pos = R['r_start'] if L['strand'] == R['strand'] or R['strand'] == '+' else R['r_end'] # Start of second segment (0-based, strand-specific)

    # Base at the break (Anchor Base)
    # For L, it is the last aligned base.
    ref_base = get_sequence(fasta_handle, bp1_chrom, bp1_pos, bp1_pos + 1) or 'N' # note this is the base of the right end of the left segment, before the event.

    # Set SV data for now. This can change based on later analysis.
    sv = {
        'chrom': bp1_chrom,        
        'pos': bp1_pos, # VCF 1-based POS is usually this value (0-based index of base before)
        'chrom2': bp2_chrom,
        'pos2': bp2_pos,
        'strands': L['strand']+R['strand'] if L['is_primary'] else R['strand']+L['strand'],
        'id': '.',
        'ref': ref_base,
        'alt': '.',
        'qual': '.',
        'filter': '.',
        'info': {}
    }
    read_seq = query_seq[read_start : read_end + 1]
    sv['info']['SEQ'] = read_seq if read_strand == '+' else reverse_complement(read_seq)

    # CASE 1: BND (Translocation OR Opposite Orientation Junction)
    # If chromosomes differ OR strands differ (opposite orientations), treat as BND.
    if bp1_chrom != bp2_chrom or L['strand'] != R['strand']:
        sv['type'] = 'BND'
        alt_seq = query_seq[L['q_end'] : R['q_start'] - 1] if read_gap else ''

        # Primary Alignment Logic for BND:
        # We must output the record anchored at the Primary Alignment.
        if L['is_primary']:

            if L['strand'] == '-':
                alt_seq = reverse_complement(alt_seq)

            # Anchor at bp1 (End of L)
            sv['alt'] = format_bnd(ref_base, alt_seq, bp2_chrom, bp2_pos, "left", L['strand'], R['strand'])
            sv['info']['SVTYPE'] = 'BND'
            sv['info']['CHR2'] = bp2_chrom
            sv['info']['POS2'] = bp2_pos

        else:
            # Anchor at bp2 (Start of R)
            # We need the anchor base at bp2_pos - 1
            if R['strand'] == '-':
                alt_seq = reverse_complement(alt_seq)

            if fasta_handle:
                anchor_R = get_sequence(fasta_handle, bp2_chrom, bp2_pos, bp2_pos + 1) or 'N'
            else:
                anchor_R = 'N' 

            sv['chrom'] = bp2_chrom
            sv['pos'] = bp2_pos
            sv['chrom2'] = bp1_chrom
            sv['pos2'] = bp1_pos
            sv['ref'] = anchor_R
            sv['alt'] = format_bnd(anchor_R, alt_seq, bp1_chrom, bp1_pos, "right", R['strand'], L['strand'])
            sv['info']['SVTYPE'] = 'BND'
            sv['info']['CHR2'] = bp1_chrom
            sv['info']['POS2'] = bp1_pos
            
        return sv

    # Continue if not a BND

    # Calculate reference span
    ref_span = bp2_pos - bp1_pos - 1
    
    # CASE 3: Duplication
    if ref_span <= 0 and read_gap >= 0:

        if bp2_pos < bp1_pos:
            sv['chrom'] = bp2_chrom
            sv['pos'] = bp2_pos - 1 # need to adjust to the base before because bp2_pos is the start of the R segment (0-based)
            sv['chrom2'] = bp1_chrom
            sv['pos2'] = bp1_pos # this is the end of the L segment, so is correct.
            ref_base = get_sequence(fasta_handle, bp2_chrom, bp2_pos - 1, bp2_pos) # This is the base before the event, 0-based.
            sv['ref'] = ref_base

        sv['type'] = 'DUP' if read_gap == 0 else 'INS'
        sv['info']['SVTYPE'] = 'DUP' if read_gap == 0 else 'INS'
        sv['info']['SVLEN'] = abs(ref_span)
        sv['info']['CHR2'] = sv['chrom2']        
        sv['info']['POS2'] = sv['pos2']        
        sv['info']['END'] = sv['pos2']        

        if fasta_handle:
            alt_seq = ref_base + query_seq[L['q_end'] : R['q_start']] if read_gap > 0 else ref_base
            # Fetch duplicated sequence
            dup_seq = get_sequence(fasta_handle, bp1_chrom, bp2_pos, bp1_pos + 1) # using first part of R segment and last of L is correct here to get the exact duplicated sequence. +1 to bp1_pos because it is 0-based and need to specify the end coordinate. 
            if dup_seq:
                sv['alt'] = alt_seq + dup_seq
            else:
                sv['alt'] = '<DUP>' if read_gap == 0 else '<INS>'
        else:
            sv['alt'] = '<DUP>' if read_gap == 0 else '<INS>'
            
        return sv

    # CASE 4: DEL
    elif ref_span > 0:

        # DELETION
        sv['type'] = 'DEL'
        sv['info']['SVTYPE'] = 'DEL'
        sv['info']['SVLEN'] = -ref_span
        sv['info']['CHR2'] = sv['chrom2']        
        sv['info']['POS2'] = sv['pos2']        
        sv['info']['END'] = sv['pos2']        

        # REF: Anchor + Deleted Sequence
        # ALT: Anchor
        if fasta_handle:
            del_seq = get_sequence(fasta_handle, bp1_chrom, bp1_pos + 1, bp2_pos) # Sequence strictly between            
            if del_seq:

                # Sanity check to see if del_seq is right length
                if len(del_seq) != abs(ref_span):
                    print("error",file=sys.stderr)

                sv['ref'] = ref_base + del_seq

            if read_gap == 0:
                sv['alt'] = ref_base

            elif read_gap > 0:
                alt_base = query_seq[L['q_end'] : L['q_end']+read_gap] 
                if L['strand'] == '-':
                    alt_base = reverse_complement(alt_base)

                sv['alt'] = ref_base + alt_base
            
            else:
                sv['alt'] = '<DEL>'
        else:
             sv['alt'] = '<DEL>'
        
        return sv

    return None

# =============================================================================
# 1. CIGAR Operations Logic
# =============================================================================

def get_cigar_indel_vcf(read, fasta_file, target_positions, target_index=None):
    """
    Detects complex indels directly from the CIGAR string.
    """
    cigar = read.cigartuples
    if not cigar: return None

    if isinstance(target_positions, int):
        target_positions = [target_positions]

    # Identify indices of all Indels (I or D)
    indel_indices = [i for i, (op, length) in enumerate(cigar) if op in [BAM_CINS, BAM_CDEL]]
    if not indel_indices: return None

    first_indel_idx = indel_indices[0]
    last_indel_idx = indel_indices[-1]
    left_flank_idx = first_indel_idx - 1
    
    # Must have a left flank that consumes reference (Match/Eq/Diff)
    if left_flank_idx < 0 or cigar[left_flank_idx][0] not in [BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF]:
        return None

    current_ref = read.reference_start
    current_read = 0
    anchor_ref_pos, anchor_read_pos = None, None
    end_ref_pos, end_read_pos = None, None

    for i, (op, length) in enumerate(cigar):
        ref_consumed = length if op in CONSUMES_REF else 0
        read_consumed = length if op in CONSUMES_READ else 0
        
        if i == left_flank_idx:
            anchor_ref_pos = current_ref + ref_consumed - 1
            anchor_read_pos = current_read + read_consumed - 1
            
        current_ref += ref_consumed
        current_read += read_consumed
        
        if i == last_indel_idx:
            end_ref_pos = current_ref
            end_read_pos = current_read
            break 

    if anchor_ref_pos is None or end_ref_pos is None: return None

    try:
        # Construct VCF record info
        out_dict = {}
        out_dict['read'] = read.query_name
        out_dict['chrom'] = read.reference_name
        out_dict['pos'] = anchor_ref_pos
        out_dict['distance'] = min(abs(out_dict['pos'] - x) for x in target_positions)
        out_dict['ref'] = fasta_file.fetch(read.reference_name, anchor_ref_pos, end_ref_pos).upper()
        out_dict['alt'] = read.query_sequence[anchor_read_pos : end_read_pos].upper()
        out_dict['chrom2'] = read.reference_name
        out_dict['pos2'] = out_dict['pos'] + len(out_dict['ref'])
        out_dict['distance2'] = min(abs(out_dict['pos2'] - x) for x in target_positions)
        out_dict['strands'] = '++'
        out_dict['type'] = 'INDEL'
        out_dict['info'] = {'Source':'CIGAR','Type':"COMPLEX" if len(indel_indices) > 1 else ("INS" if cigar[indel_indices[0]][0] == BAM_CINS else "DEL"),'Read':read.query_name,'OrigCigar':read.cigarstring,'ReadSeq':read.query_sequence}
        return out_dict
    except (ValueError, IndexError):
        return None

# =============================================================================
# 2. Supplementary Alignment (SA) Logic
# =============================================================================

def get_sa_indel_vcf(read, fasta_file, target_positions, target_index):
    """
    Detects events by merging Primary and Supplementary Alignments.
    Handles DEL/INS for collinear events and BND for translocations/inversions.
    """
    if not read.has_tag("SA"): return None

    sa_tag = read.get_tag('SA')
        
    # 1. Get Read Seq
    query_seq = read.query_sequence
    if read.is_reverse:
        query_seq = reverse_complement(query_seq)
        
    # 2. Parse Primary Alignment (ALN1) - Mark as PRIMARY
    cigar1 = read.cigar
    aln1 = cigar_to_coords(
        cigar1, 
        read.reference_name, 
        read.reference_start, # start position of alignment, 0-based
        '-' if read.is_reverse else '+',
        is_primary=True
    )

    # 3. Parse Secondary Alignment(s) from SA tag (ALN2)
    sa_records = [s for s in sa_tag.split(';') if s]
    
    if not sa_records: return None
    
    sa = sa_records[0]
    parts = sa.split(',')
    if len(parts) < 4:
        print(f"Error parsing SA tag {sa} in read {read.query_name}", file=sys.stderr)
        return None
    
    rname_2 = parts[0]
    pos_2 = int(parts[1]) - 1 # pos is 1-based in SA tags. Convert to 0-based. This is the left-most mapping pos.
    strand_2 = parts[2]
    cigar_str_2 = parts[3]
    sa_mapq = int(parts[4])
    sa_nm = int(parts[5])
    
    cigar2 = parse_cigar_string(cigar_str_2)
    
    aln2 = cigar_to_coords(
        cigar2,
        rname_2,
        pos_2,
        strand_2,
        is_primary=False
    )
        
    # The cigars are not aligned if aln1 and aln2 are different directions
    # so flip the aln2 coordinates on the read.
    if aln1['strand'] != aln2['strand']:
        if aln1['strand'] == '-':
            q_end = read.query_length - aln1['q_start'] - 1
            q_start = read.query_length - aln1['q_end'] - 1
            aln1['q_start'] = q_start
            aln1['q_end'] = q_end
        else:
            q_end = read.query_length - aln2['q_start'] - 1
            q_start = read.query_length - aln2['q_end'] - 1
            aln2['q_start'] = q_start
            aln2['q_end'] = q_end
    
    # 4. Call SV
    sv = call_sv_from_split_read(aln1, aln2, query_seq, fasta_handle=fasta_file)
    
    if sv:
        # Add detailed read info to tagsx
        sv['read'] = read.query_name
        sv['info']['RN'] = read.query_name
        sv['info']['STRAND'] = aln1['strand']
        sv['info']['CIGAR'] = read.cigarstring
        # Sanitize SA tag (replace semicolons to preserve VCF format)
        sv['info']['SA'] = sa.replace(';', '')
        sv['info']['SAMAPQ'] = sa_mapq

        sv['distance'] = min(abs(sv['pos'] - x) for x in target_positions)
        sv['distance2'] = target_index.get_min_dist(sv['chrom2'],sv['pos2'])# min(abs(sv['pos2'] - x) for x in target_positions)
                        
        # Common Info Fields
        sv['info']['Source'] = 'SA'
        
        return sv

    return None

# =============================================================================
# 3. Soft Clip Logic
# =============================================================================

def get_softclip_indel_vcf(read, fasta_file, target_positions, search_range, min_clip=8):
    """
    Detects indels by realigning soft-clipped sequences to the reference.
    """
    if read.is_supplementary or read.has_tag("SA") or not read.cigartuples:
        return None

    # get primary alignment
    cigar1 = read.cigar
    aln1 = cigar_to_coords(
        cigar1, 
        read.reference_name, 
        read.reference_start, # start position of alignment, 0-based
        '-' if read.is_reverse else '+',
        is_primary=True
    )

    clip_side = None
    clip_seq_str = None
    clip_len = 0
    ref_start = 0
    ref_end = 0
    clip_read_start = 0
    clip_read_end = 0

    # Check Right Clip
    if cigar1[-1][0] == BAM_CSOFT_CLIP and cigar1[-1][1] >= min_clip:
        clip_side = 'Right'
        clip_len = cigar1[-1][1]
        clip_seq_str = read.query_sequence[-clip_len:]
        ref_start = read.reference_end
        clip_read_start = read.query_alignment_end

    # Check Left Clip
    elif cigar1[0][0] == BAM_CSOFT_CLIP and cigar1[0][1] >= min_clip:
        clip_side = 'Left'
        clip_len = cigar1[0][1]
        clip_seq_str = read.query_sequence[:clip_len]
        ref_end = read.reference_start
        clip_read_end = read.query_alignment_start

    else:
        return None

    # Get search window origin
    clip_ref_pos = read.reference_end if clip_side == 'Right' else read.reference_start
    window_start_genomic = clip_ref_pos - search_range
    
    query_chunk = clip_seq_str
    try:
        ref_window = fasta_file.fetch(read.reference_name, window_start_genomic, clip_ref_pos + search_range).upper()
    except ValueError:
        return None

    # Simple exact match search first
    aln2 = None
    match_count = ref_window.count(query_chunk)
    match_index = ref_window.find(query_chunk)
    
    if match_count > 1:
        return None
    
    elif match_count == 1 and match_index > 0:
        if clip_side == 'Right':
            aln2 = cigar_to_coords([(4,read.query_length - clip_len),(0,len(query_chunk))], 
                                   read.reference_name, 
                                   window_start_genomic + match_index, 
                                   '-' if read.is_reverse else '+',
                                   False)
        else: # Left
            aln2 = cigar_to_coords([(0,len(query_chunk)),(4,read.query_length - clip_len)], 
                                   read.reference_name, 
                                   window_start_genomic + match_index, 
                                   '-' if read.is_reverse else '+',
                                   False)
        
    else:
        result = edlib.align(query_chunk, ref_window, mode="HW", task="path", k=-1)
        if result['locations'] and len(result['locations']) == 0:
            match_cigar = parse_cigar_string(result['cigar'])
            if clip_side == "Right" and match_cigar[-1][0] == BAM_CMATCH and match_cigar[-1][1] >= min_clip: # get ref coordinate of right-most anchor of the match 
                aln2 = cigar_to_coords([(4,read.query_length - match_cigar[-1][1]),(0,match_cigar[-1][1])], 
                                       read.reference_name, 
                                       window_start_genomic + result['locations'][0][1] - match_cigar[-1][1] + 1,
                                       '-' if read.is_reverse else '+',
                                       False)

            elif clip_side == "Left" and match_cigar[0][0] == BAM_CMATCH and match_cigar[0][1] >= min_clip: # get ref coordinate of left
                aln2 = cigar_to_coords([(0,match_cigar[-1][1]),(4,read.query_length - match_cigar[0][1])],
                                         read.reference_name,
                                         window_start_genomic + result['locations'][0][0],
                                         '-' if read.is_reverse else '+',
                                         False)
            
            else:
                return None
        else:
            return None

    query_seq = read.query_sequence
    if read.is_reverse:
        query_seq = reverse_complement(query_seq)

    sv = call_sv_from_split_read(aln1, aln2, query_seq, fasta_handle=fasta_file)

    if sv:
        sv['read'] = read.query_name
        sv['distance'] = min(abs(sv['pos'] + 1 - x) for x in target_positions)
        sv['distance2'] = min(abs(sv['pos2'] - x) for x in target_positions)
        sv['info']['Source'] = 'SoftClip'
    
        return sv

    return None


# ============================================================================
# SECTION 1: ORIGINAL INDEL CALCULATION FUNCTIONS
# ============================================================================

# function to count the number of reads that support an indel or BND event
def add_normal_counts(df, reads, fasta, handicap=5,window=25,flank=300):

    # make ref and alt sequences for each indel/BND
    results = []

    for _, row in df.iterrows():
        # 0-based index of the variant start
        # VCF POS is 1-based, so pos-1 is the 0-based index of the anchor base
        start_idx = int(row['pos']) - 1
        
        # 1. Get REF sequence with flanks
        # Interval: [pos - 1 - flank, pos - 1 + len(ref) + flank]
        # This covers the anchor, the rest of the ref allele, and both flanks.
        ref_seq = fasta.fetch(
            row['chrom'], 
            start_idx - flank, 
            start_idx + len(row['ref']) + flank
        )
        
        alt_seq = ref_seq
        if row['alt'] != '.':
            alt_seq = generate_contig(row['chrom'], row['pos'], row['ref'], row['alt'], row['type'], fasta, flank)
        
        results.append({
            'refseq': ref_seq,
            'altseq': alt_seq,
            'control_alt_counts': 0,
            'control_total_counts': 0
        })

    # Assign back to DataFrame efficiently
    result_df = pd.DataFrame(results)
    df = pd.concat([df.reset_index(drop=True), result_df], axis=1)

    total_reads = set()

    # iterate through reads in bam file for the first position
    for read in reads:
        # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
        if read.is_mapped is False or \
            read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.mapping_quality == 0:
            continue

        total_reads.add(read.query_name)

        # for speed: if the read has no indels, add it to the ref_count set
        if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == len(read.query_sequence) and read.has_tag('SA') is False:
            continue
        
        for it, row in df.iterrows():                        

            # if cigar has any D/I operations, get indel from read
            if ('D' in read.cigarstring or 'I' in read.cigarstring):

                vcf_dict = get_cigar_indel_vcf(read, fasta, row['pos'])

                if vcf_dict and vcf_dict['pos'] == row['pos'] and vcf_dict['ref'] == row['ref'] and vcf_dict['alt'] == row['alt']:
                    df.at[it,'control_alt_counts'] += 1
                    continue

            read_seq = read.query_sequence
            
            if read_seq in row['refseq']:
                continue

            if read_seq in row['altseq']:
                df.at[it,'control_alt_counts'] += 1
                continue

            ref_align = edlib.align(read_seq, row['refseq'], mode="HW", task="path") #align.align_optimal(seq.NucleotideSequence(row['refseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]
            alt_align = edlib.align(read_seq, row['altseq'], mode="HW", task="path") #align.align_optimal(seq.NucleotideSequence(row['altseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]

            if alt_align['editDistance'] < ref_align['editDistance']: #.score > 0 and alt_align.score - handicap > ref_align.score:
                ref_align_cigar_sum = cigar_summary(ref_align['cigar'])
                alt_align_cigar_sum = cigar_summary(alt_align['cigar'])

                if ((len(row['ref']) > len(row['alt']) and ref_align_cigar_sum['D'] >= alt_align_cigar_sum['D']) or
                    (len(row['ref']) < len(row['alt']) and ref_align_cigar_sum['I'] >= alt_align_cigar_sum['I']) or
                    len(row['ref']) == len(row['alt'])):
                    df.at[it,'control_alt_counts'] += 1
                    # print to stderr: found control alt count
                    print(f"\tFound control alt count for {row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{read.query_name}", file=sys.stderr)

    df['control_total_counts'] = len(total_reads)
    return df.copy()

# ============================================================================
# SECTION 2: CRISPR PREDICTION FUNCTIONS
# ============================================================================

def detect_deletion_from_softclips(read, min_softclip_size=3, max_gap=50):
    """Detect potential large deletions from soft-clips and supplementary alignments."""
    if not read.cigartuples:
        return 0
    
    total_deletion = 0
    
    # Parse CIGAR operations
    cigar_ops = []
    for op, length in read.cigartuples:
        cigar_ops.append((op, length))
    
    # Count soft-clips as potential deletions
    for op, length in cigar_ops:
        if op == 4:  # Soft-clip
            if length >= min_softclip_size:
                total_deletion += length
    
    # Look for patterns: softclip -> deletion -> softclip
    for i in range(len(cigar_ops) - 2):
        op1, len1 = cigar_ops[i]
        op2, len2 = cigar_ops[i + 1]
        op3, len3 = cigar_ops[i + 2]
        
        if op1 == 4 and op2 == 2 and op3 == 4:  # Softclip -> Deletion -> Softclip
            if len1 >= min_softclip_size and len2 <= max_gap and len3 >= min_softclip_size:
                total_deletion += len2
        elif op1 == 2 and op2 == 4 and op3 == 2:  # Deletion -> Softclip -> Deletion
            if len2 >= min_softclip_size and len1 <= max_gap and len3 <= max_gap:
                total_deletion += (len1 + len3)
    
    # Look for single softclip patterns with nearby deletions
    for i in range(len(cigar_ops) - 1):
        op1, len1 = cigar_ops[i]
        op2, len2 = cigar_ops[i + 1]
        
        if op1 == 4 and op2 == 2:  # Softclip -> Deletion
            if len1 >= min_softclip_size and len2 <= max_gap:
                total_deletion += len2
        elif op1 == 2 and op2 == 4:  # Deletion -> Softclip
            if len2 >= min_softclip_size and len1 <= max_gap:
                total_deletion += len1
    
    # Look for end soft-clips (left/right) representing large deletions
    if len(cigar_ops) >= 1:
        # Left soft-clip
        if cigar_ops[0][0] == 4:
            left_softclip = cigar_ops[0][1]
            if left_softclip >= min_softclip_size:
                middle_align = sum(length for op, length in cigar_ops[1:] if op == 0)
                if middle_align > 5:
                    total_deletion += left_softclip
        
        # Right soft-clip
        if cigar_ops[-1][0] == 4:
            right_softclip = cigar_ops[-1][1]
            if right_softclip >= min_softclip_size:
                middle_align = sum(length for op, length in cigar_ops[:-1] if op == 0)
                if middle_align > 5:
                    total_deletion += right_softclip
    
    # Look for internal soft-clips
    for i, (op, length) in enumerate(cigar_ops):
        if op == 4 and length >= min_softclip_size:
            left_align = sum(cigar_ops[j][1] for j in range(i) if cigar_ops[j][0] == 0)
            right_align = sum(cigar_ops[j][1] for j in range(i+1, len(cigar_ops)) if cigar_ops[j][0] == 0)
            
            if left_align > 2 and right_align > 2:
                total_deletion += length
    
    # Check supplementary alignments (SA tags) for soft-clip deletions
    if read.has_tag('SA'):
        sa_tag = read.get_tag('SA')
        sa_entries = sa_tag.split(';')
        
        for sa_entry in sa_entries:
            if not sa_entry:
                continue
            
            sa_parts = sa_entry.split(',')
            if len(sa_parts) >= 4:
                sa_cigar = sa_parts[3]
                sa_softclips = parse_cigar_for_softclips(sa_cigar)
                for softclip_size in sa_softclips:
                    if softclip_size >= min_softclip_size:
                        total_deletion += softclip_size
    
    return total_deletion

def parse_cigar_for_softclips(cigar_string):
    """Parse CIGAR string to extract soft-clip sizes."""
    softclips = []
    import re
    
    cigar_parts = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    
    for length_str, operation in cigar_parts:
        if operation == 'S':  # Soft-clip
            softclips.append(int(length_str))
    
    return softclips

def parse_cigar_once(read):
    """Parse CIGAR operations to extract insertion, deletion, and soft-clip counts."""
    cigar_summary = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
    if not hasattr(read, 'cigartuples') or read.cigartuples is None:
        return {'insertions': 0, 'deletions': 0, 'softclips': 0}
    
    for op, length in read.cigartuples:
        if op == 1:
            cigar_summary['I'] += length
        elif op == 2:
            cigar_summary['D'] += length
        elif op == 4:
            cigar_summary['S'] += length
    
    # Total deletions = explicit deletions + deletions from soft-clips
    total_deletions = cigar_summary['D'] + detect_deletion_from_softclips(read)
    
    return {
        'insertions': cigar_summary['I'],
        'deletions': total_deletions,
        'softclips': cigar_summary['S']
    }

def count_mismatches_fast(read):
    """
    Returns the number of single-base mismatches (SNPs) in the alignment.
    Excludes Indels and Soft-clips.
    """
    try:
        # 1. Get total Edit Distance (Mismatches + Indels)
        nm_tag = read.get_tag("NM")
        
        # 2. Get CIGAR statistics
        # read.get_cigar_stats() returns a tuple. The first element is a list of counts 
        # by CIGAR operation (M, I, D, N, S, H, P, =, X, B).
        # Index 1 = Insertions (I)
        # Index 2 = Deletions (D)
        stats = read.get_cigar_stats()[0]
        insertion_len = stats[1]
        deletion_len = stats[2]
        
        # 3. Subtract Indels from NM to isolate Mismatches
        return nm_tag - (insertion_len + deletion_len)
        
    except (KeyError, IndexError):
        # Handle reads without NM tags or unmapped reads
        return 0

def calculate_control_fractions(control_bam, chrom, position, window=50):
    """Calculate fractions of control reads with different variant types."""
    fractions = {
        'fraction_control_reads_del': 0.0,
        'fraction_control_reads_ins': 0.0,
        'fraction_control_reads_mismatch': 0.0,
        'fraction_control_reads_softclip': 0.0
    }
    if not control_bam:
        return fractions
        
    total = del_count = ins_count = mismatch_count = softclip_count = 0
    
    for read in control_bam.fetch(chrom, position - window, position + window):
        if read.is_unmapped or read.is_duplicate:
            continue
        total += 1
        if read.cigartuples:
            if any(op == 2 for op, _ in read.cigartuples): 
                del_count += 1
            if any(op == 1 for op, _ in read.cigartuples): 
                ins_count += 1
            if any(op == 4 for op, _ in read.cigartuples): 
                softclip_count += 1
        if read.has_tag('MD') and any(c.isalpha() for c in read.get_tag('MD')):
            mismatch_count += 1
    
    if total > 0:
        fractions.update({
            'fraction_control_reads_del': del_count / total,
            'fraction_control_reads_ins': ins_count / total,
            'fraction_control_reads_mismatch': mismatch_count / total,
            'fraction_control_reads_softclip': softclip_count / total
        })
    return fractions

def calculate_exclusivity_features(read, control_bam, chrom, position, window=100):
    """Calculate whether variants are exclusive to edited samples."""
    exclusivity = {
        'variant_exclusive_to_edited': 0,
        'insertion_exclusive_to_edited': 0,
        'deletion_exclusive_to_edited': 0,
        'control_has_same_variant': 0,
    }
    if not control_bam or read.is_unmapped or read.is_duplicate or not read.cigartuples:
        return exclusivity

    has_insertion = any(op == 1 for op, _ in read.cigartuples)
    has_deletion = any(op == 2 for op, _ in read.cigartuples)
    
    read_start, read_end = read.reference_start, read.reference_end

    control_has_same_ins = False
    control_has_same_del = False
    
    for control_read in control_bam.fetch(chrom, position - window, position + window):
        if control_read.is_unmapped or control_read.is_duplicate or not control_read.cigartuples:
            continue
        ctrl_start, ctrl_end = control_read.reference_start, control_read.reference_end
        
        # Check for overlapping reads
        if not (read_end < ctrl_start or ctrl_end < read_start):
            if has_insertion and any(op == 1 for op, _ in control_read.cigartuples):
                control_has_same_ins = True
            if has_deletion and any(op == 2 for op, _ in control_read.cigartuples):
                control_has_same_del = True

    insertion_exclusive = has_insertion and not control_has_same_ins
    deletion_exclusive = has_deletion and not control_has_same_del
    
    any_exclusive = insertion_exclusive or deletion_exclusive
    control_has_similar = control_has_same_ins or control_has_same_del

    exclusivity['deletion_exclusive_to_edited'] = 1 if deletion_exclusive else 0
    exclusivity['control_has_same_variant'] = 1 if control_has_similar else 0
    return exclusivity

def calculate_distance_to_closest_pam(position, targets_df, chrom):
    """Calculate distance to closest PAM site."""
    if targets_df is None or len(targets_df) == 0:
        return -1
    distances = np.abs(targets_df['Start'] - position)
    return distances.min() if len(distances) > 0 else -1

def predict_reads_at_position(bam_file, chrom, start, end, pampos, model, fasta, is_on_target=0, control_bam=None, threshold=0.80):
    
    reads = []
    for read in bam_file.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate:
            continue
        reads.append(read)
    if not reads:
        return 0, 0.0

    # Extract features for each read
    features_list = []
    for read in reads:
        cigar_data = parse_cigar_once(read)
        
        read_start = read.reference_start
        read_end = read.reference_end
        
        # Check if read overlaps with target sites
        is_at_any_target = 0
        if pampos is not None:
            for p in pampos:
                if read_start <= p + 25 and read_end >= p - 25:
                    is_at_any_target = 1

        
        # Get control fractions
        control_fractions = calculate_control_fractions(control_bam, chrom, start, window=50) if control_bam else {
            'fraction_control_reads_del': 0.0,
            'fraction_control_reads_ins': 0.0,
            'fraction_control_reads_mismatch': 0.0,
            'fraction_control_reads_softclip': 0.0
        }
        
        # Calculate exclusivity features
        exclusivity = calculate_exclusivity_features(read, control_bam, chrom, start, window=100) if control_bam else {
            'deletion_exclusive_to_edited': 0,
            'control_has_same_variant': 0,
        }
        
        # Convert read features to binary
        read_has_deletion = 1 if cigar_data['deletions'] > 0 else 0
        read_has_insertion = 1 if cigar_data['insertions'] > 0 else 0
        read_has_mismatch = 1 if count_mismatches_fast(read) > 0 else 0
        
        # Calculate indel characteristics
        total_indel_size = cigar_data['insertions'] + cigar_data['deletions']
        
        # Indel size category
        indel_size_category = 0
        if total_indel_size > 0:
            if total_indel_size <= 3:
                indel_size_category = 1
            elif total_indel_size <= 10:
                indel_size_category = 2
            else:
                indel_size_category = 3
        
        # Insertion to deletion ratio
        insertion_to_deletion_ratio = 0.0
        if cigar_data['deletions'] > 0:
            insertion_to_deletion_ratio = cigar_data['insertions'] / cigar_data['deletions']
        elif cigar_data['insertions'] > 0:
            insertion_to_deletion_ratio = 10.0
        else:
            insertion_to_deletion_ratio = 0.0
        
        # Indel complexity score
        indel_complexity_score = 0.0
        if read.cigartuples:
            indel_operations = [op for op, length in read.cigartuples if op in [1, 2]]
            complexity = len(indel_operations) + (total_indel_size / 10.0)
            indel_complexity_score = min(complexity, 10.0)
        
        # Create features dictionary
        features = {
            'read_pair_gap': abs(read.template_length) if hasattr(read, 'template_length') and read.is_paired and read.is_proper_pair else -1,
            'read_insertion': cigar_data['insertions'],
            'read_deletion': cigar_data['deletions'],
            'read_mismatch': count_mismatches_fast(read),
            'read_softclip': cigar_data['softclips'],
            'read_del_vs_control': read_has_deletion - control_fractions['fraction_control_reads_del'],
            'read_ins_vs_control': read_has_insertion - control_fractions['fraction_control_reads_ins'],
            'read_mismatch_vs_control': read_has_mismatch - control_fractions['fraction_control_reads_mismatch'],
            'deletion_exclusive_to_edited': exclusivity['deletion_exclusive_to_edited'],
            'control_has_same_variant': exclusivity['control_has_same_variant'],
            'distance_to_closest_pam': min(abs(x-start) for x in pampos) if pampos is not None else -1,
            'is_on_target_site': is_on_target,
            'is_at_any_target_site': is_at_any_target,
            'total_indel_size': total_indel_size,
            'indel_size_category': indel_size_category,
            'insertion_to_deletion_ratio': insertion_to_deletion_ratio,
            'indel_complexity_score': indel_complexity_score
        }
        
        features_list.append(features)
    
    if not features_list:
        return 0, 0.0

    # Get expected features for the model
    if hasattr(model, 'feature_names_in_'):
        expected_features = list(model.feature_names_in_)
    else:
        expected_features = ['read_pair_gap', 'read_insertion', 'read_deletion', 'read_mismatch', 'read_softclip', 'read_del_vs_control', 'read_ins_vs_control', 'read_mismatch_vs_control', 'deletion_exclusive_to_edited', 'control_has_same_variant', 'is_on_target_site', 'is_at_any_target_site', 'distance_to_closest_pam', 'total_indel_size', 'indel_size_category', 'insertion_to_deletion_ratio', 'indel_complexity_score']
    
    features_df = pd.DataFrame(features_list)
    features_df = features_df.reindex(columns=expected_features, fill_value=0)
    
    # Get model predictions
    preds = model.predict_proba(features_df)[:, 1]
    
    # Calculate results
    avg_probability = float(preds.mean() * 100)
    
    return int((preds >= threshold).sum()), avg_probability

def is_in_target(chrom, pos_0based, lookup):
    if chrom not in lookup: return False
    starts, ends = lookup[chrom]
    
    # Find index where pos would insert into starts
    idx = bisect_right(starts, pos_0based)
    
    if idx == 0: return False
    
    # Check if the interval at idx-1 actually covers the position
    # (This assumes non-overlapping intervals)
    return ends[idx-1] > pos_0based

# ============================================================================
# SECTION 3: MAIN FUNCTION
# ============================================================================

def main():
    
    parser = argparse.ArgumentParser(description='Extract variant reads and predict CRISPR reads')

    # BAM files, target VCF, and reference FASTA.
    parser.add_argument('--edited-bam',type=str,required=True,help='Edited/experimental BAM/CRAM file')
    parser.add_argument('--control-bam',type=str,required=True,help='Control BAM/CRAM file')
    parser.add_argument('--target-file',type=str,required=True,help='Target sites coordinate file (CSV/BED format) - used for both bed coordinates and CRISPR targets')
    parser.add_argument('-f','--fasta',type=str,default="/storage2/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh_v4.3.6/hg38_PLVM_CD19_CARv4_cd34.fa",help='Reference fasta file')

    # Search and filtering parameters
    parser.add_argument('-w','--target-window',type=int,default=150,help='Window size')
    parser.add_argument('-d','--max-mutation-distance',type=int,default=25,help='Distance')
    parser.add_argument('-s','--mutation-search-window',type=int,default=20000,help='Search window')
    parser.add_argument('-l','--min-softclip-length',type=int,default=8,help='Minimum softclip length')
    parser.add_argument('-b','--min-bnd-mapqual',type=int,default=40,help='Minimum BND mapping quality')
    parser.add_argument('-m','--min-coverage',type=int,default=1,help='Minimum reads')
    parser.add_argument('-x','--max-in-control',type=int,default=5,help='Maximum supporting reads in control/unedited sample to report an indel/bnd event.')
    parser.add_argument('-q','--min-mapqual',type=int,default=20,help='Minimum mapping quality')
    parser.add_argument('-n','--max-read-mismatches',type=int,default=4,help='Maximum number of mismatches')

    # Optionally search at targets from one chromosome
    parser.add_argument('-c','--chromosome',type=str,default=None,help='Chromosome to process')

    # Outputs
    parser.add_argument('-o','--outfile',type=str,help='Output file (optional)')
    parser.add_argument('-u','--unevaluable-reads-logfile',type=str,help='File with information on reads that were not evaluable.')
    parser.add_argument('-v','--verbose',action='store_true',help='Print verbose output')
    
    # CRISPR prediction arguments
    parser.add_argument('--crispr-model',type=str,default='models/site14_site5_combined_model.pkl',help='Trained CRISPR ML model file (.pkl) for read prediction (default: site14_site5_combined_model.pkl)')
    parser.add_argument('--crispr-threshold',type=float,default=0.70,help='Probability threshold for CRISPR prediction (default: 0.70)')
    parser.add_argument('--enable-crispr-prediction',action='store_true',help='Enable CRISPR read prediction (uses default model and target-file as targets)')
    parser.add_argument('--targets-csv',type=str,help='Alternative target sites CSV file for feature extraction (optional - uses target-file if not specified)')
    parser.add_argument('--filter-off-target-fp', action='store_true', help='Filter off-target sites that are likely false positives (indel reads > 0 but no predicted CRISPR reads).')
    parser.add_argument('--fp-log', type=str, help='Log file for filtered false positive off-target sites.')
    
    args = parser.parse_args()

    if args.verbose:
        print("Processing input file", file=sys.stderr)

    # ========================================================================
    # STEP 1: Process input BED file and create genomic intervals
    # ========================================================================

    # target regions, in VCF format.
    vcf_data = []

    try:
        # Open VCF file with pysam
        vcf_in = pysam.VariantFile(args.target_file)
        
        for rec in vcf_in:
            start = rec.pos - 1
            end = rec.pos
            info_dict = dict(rec.info)
            is_target = 1 if info_dict.get('TARGET', False) else 0
            
            vcf_data.append({
                'Chromosome': rec.chrom,
                'Start': start,
                'End': end,
                'Pos': rec.pos,  # Keep 1-based POS for reference/output
                'Info': info_dict,
                'Ontarget': is_target
            })

    except Exception as e:
        print(f"Error reading VCF input: {e}", file=sys.stderr)

    # Create DataFrame from list of dicts
    bedDf = pd.DataFrame(vcf_data)

    # Generate coordinate lookup index:
    target_index = GenomicDistanceIndex(bedDf)

    # Check for empty input data
    if len(bedDf) == 0:
        print("Warning: Input target file has no data rows. Writing empty output.", file=sys.stderr)
        # Write empty output file with header only
        output_columns = ['Cluster', 'Chromosome', 'Start', 'End', 'Ontarget', 'Gene', 'indel_type', 
                        'indel_fraction', 'indel_allele_fraction', 'indel_size', 'indel_bases', 
                        'num_edited', 'num_control', 'total_edited', 'total_control', 'significance',
                        'prediction', 'probability', 'model_info']
        empty_df = pd.DataFrame(columns=output_columns)
        if args.outfile:
            empty_df.to_csv(args.outfile, sep='\t', index=False)
        else:
            empty_df.to_csv(sys.stdout, sep='\t', index=False)
        sys.exit(0)

    # Create PyRanges object and cluster intervals
    bedPr = pr.PyRanges(bedDf[['Chromosome', 'Start', 'End', 'Pos', 'Info', 'Ontarget']])
    bedPr = bedPr.cluster(slack=args.target_window)
    mergedBedPr = bedPr.merge(by='Cluster', strand=False, slack=args.target_window)

    # Join aggregated data back to the merged intervals
    mergedBedDf = mergedBedPr.df.join(
        bedPr.df.groupby('Cluster')['Pos'].agg(list).reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )
    mergedBedDf = mergedBedDf.join(
        bedPr.df.groupby('Cluster')['Info'].agg(list).reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )
    mergedBedDf = mergedBedDf.join(
        bedPr.df.groupby('Cluster')['Ontarget'].agg('max').reset_index().set_index('Cluster'),
        on='Cluster', how='left'
    )

    # Filter by chromosome if specified
    if args.chromosome is not None:
        print("Processing chromosome", args.chromosome, file=sys.stderr)
        if args.chromosome.startswith('chr'):
            mergedBedDf = mergedBedDf[mergedBedDf['Chromosome'] == args.chromosome]
        else:
            mergedBedDf = mergedBedDf[(mergedBedDf['Chromosome'] == args.chromosome) | 
                                    (mergedBedDf['Chromosome'] == f"chr{args.chromosome}")]
        
        if len(mergedBedDf) == 0:
            print(f"Warning: No data found for chromosome {args.chromosome}", file=sys.stderr)
            all_chromosomes = sorted(bedPr.df['Chromosome'].unique())
            print(f"Available chromosomes: {all_chromosomes}", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"Found {len(mergedBedDf)} intervals for chromosome {args.chromosome}", file=sys.stderr)

    if args.verbose:
        print("Done processing input file", file=sys.stderr)

    # ========================================================================
    # STEP 2: Open BAM files and reference
    # ========================================================================
    
    edited_bamfile = pysam.AlignmentFile(args.edited_bam,"rc",reference_filename=args.fasta)
    control_bamfile = pysam.AlignmentFile(args.control_bam,"rc",reference_filename=args.fasta)
    refFasta = pysam.FastaFile(args.fasta)

    # make outfile to print to, or use stdout
    if args.outfile:
        fp = open(args.outfile, 'w')
    else:
        fp = sys.stdout

    # ========================================================================
    # STEP 3: Load CRISPR prediction model and targets (if enabled)
    # ========================================================================    

    crispr_model = None
    if args.enable_crispr_prediction:
        try:
            print(f"Loading CRISPR ML model: {args.crispr_model}", file=sys.stderr)
            crispr_model = joblib.load(args.crispr_model)
            print(f"CRISPR model loaded successfully", file=sys.stderr)
            print(f"CRISPR prediction enabled with threshold {args.crispr_threshold}", file=sys.stderr)
        except Exception as e:
            print(f"Error loading CRISPR model '{args.crispr_model}': {e}", file=sys.stderr)
            print(f"Make sure the model file exists in the current directory or provide full path", file=sys.stderr)
            sys.exit(1)

    # ========================================================================
    # STEP 4: Create output header
    # ========================================================================
    
    header_columns = ['chrom', 'start', 'end', 'pam_positions', 'total_reads', 'indel_reads', 'indel_fraction', 
                     'control_reads', 'control_indel_reads', 'control_indel_fraction', 'indel_count', 'indel_info', 
                     'bnd_count', 'bnd_info', 'target_info', 'is_target']
    if args.enable_crispr_prediction:
        header_columns.extend(['crispr_predicted_reads', 'crispr_prediction_fraction', 'crispr_prediction_probability'])

    print("\t".join(header_columns), file=fp, flush=True)

    fp_log = None
    if args.fp_log:
        fp_log = open(args.fp_log, 'w')
        fp_log.write("\t".join(header_columns) + "\n")

    unevaluable_read_log = None
    if args.unevaluable_reads_logfile:
        unevaluable_read_log = open(args.unevaluable_reads_logfile, 'w')

    # ========================================================================
    # STEP 5: Process each genomic interval
    # ========================================================================
    
    for _, row in mergedBedDf.iterrows():
        if args.verbose:
            print(f"Processing interval {row['Chromosome']}:{row['Start']}-{row['End']}", file=sys.stderr)

        # Set up edit df
        readaln = pd.DataFrame(columns=['read','chrom','pos','chrom2','pos2','ref','alt','strands','type','info'])
    
        if args.verbose:
            print(f"\tGetting reads that align within window of {row['Chromosome']}:{row['Start']}-{row['End']}", file=sys.stderr)

        # get reads that align within a defined region containing the merged target interval
        for read in edited_bamfile.fetch(row['Chromosome'], row['Start']-args.target_window, row['End']+args.target_window, multiple_iterators = True):

            # skip if not primary alignment or a duplicate or poor mapping quality
            if (read.is_mapped is False or
                read.is_duplicate is True or
                read.is_secondary is True or
                read.is_supplementary is True or
                read.mapping_quality < args.min_mapqual or 
                count_mismatches_fast(read) > args.max_read_mismatches):
                continue

            # Determine whether the read pair has the proper orientation, that is: ---> <---
            proper_paired_read = (read.is_paired and 
                                  read.is_forward != read.is_reverse and 
                                  ((read.is_forward and read.reference_start <= read.next_reference_start) or
                                   (read.is_reverse and read.reference_start >= read.next_reference_start)))
                                   
            cigar = read.cigartuples # get cigar info
            mate_cigar = parse_cigar_string(read.get_tag('MC')) if read.has_tag('MC') else None # get mate_cigar

            read_strand = "+"
            if read.is_reverse is True:
                read_strand = "-"

            # dict to store mutation info in VCF record format
            vcf_dict = None

            # if cigar has any D/I operations
            if ('D' in read.cigarstring or 'I' in read.cigarstring) and proper_paired_read:

                if args.verbose:
                    print(f"\tAnalyzing cigars in {read.query_name} from {row['Chromosome']}:{row['Start']-args.target_window}-{row['End']+args.target_window}", file=sys.stderr)

                vcf_dict = get_cigar_indel_vcf(read, refFasta, row['Pos'], target_index)

                # if no indel is found or its too far away from the closest PAM position
                if (vcf_dict is None or 
                    vcf_dict['distance'] > args.max_mutation_distance and vcf_dict['distance2'] > args.max_mutation_distance):
                    # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                    if unevaluable_read_log:
                        print(f"{str(read)}\t{vcf_dict}", file=unevaluable_read_log)

                    continue

            # read has supplementary alignments
            elif read.has_tag('SA'):            

                if args.verbose:
                    print(f"\tAnalyzing SA in {read.query_name}, {read.get_tag('SA') if read.has_tag('SA') else 'None'} from {row['Chromosome']}:{row['Start']-args.target_window}-{row['End']+args.target_window}", file=sys.stderr)

                vcf_dict = get_sa_indel_vcf(read, refFasta, row['Pos'], target_index)

                # skip if indel is too far away from the PAM position
                if (vcf_dict is None or 
                    (vcf_dict['type']=='INDEL' and 
                        vcf_dict['distance'] > args.max_mutation_distance and vcf_dict['distance2'] > args.max_mutation_distance)):

                    # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                    if unevaluable_read_log:
                        print(f"{str(read)}\t{vcf_dict}", file=unevaluable_read_log)

                    continue
                
                # if SA is an indel then it should be bounded by the read pair ends. If not then continue.
                if (vcf_dict['type']=='INDEL' and proper_paired_read and
                    (min([vcf_dict['pos'],vcf_dict['pos2']]) < min([read.reference_start,read.next_reference_start]) or
                     max([vcf_dict['pos'],vcf_dict['pos2']]) > min([read.reference_start,read.next_reference_start])+read.template_length)):
                        
                        # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                        if unevaluable_read_log:
                            print(f"{str(read)}\t{vcf_dict}", file=unevaluable_read_log)

                        continue
                
                # if SA is a BND, check to see if the other end is in the target list
                if (vcf_dict['type']=='BND' and 
                    (read.mapping_quality < args.min_bnd_mapqual or vcf_dict['info']['SAMAPQ'] < args.min_bnd_mapqual and 
                     vcf_dict['distance'] > args.max_mutation_distance and vcf_dict['distance2'] > args.max_mutation_distance)): 
        
                        # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                        if unevaluable_read_log:
                            print(f"{str(read)}\t{vcf_dict}", file=unevaluable_read_log)

                        continue

            # read has softclips
            elif cigar[0][0] >= 4 and cigar[0][1] >= args.min_softclip_length or \
                    cigar[-1][0] >= 4 and cigar[-1][1] >= args.min_softclip_length:

                if args.verbose:
                    print(f"\tAnalyzing softclips in {read.query_name}, {read.cigarstring}, from {row['Chromosome']}:{row['Start']-args.target_window}-{row['End']+args.target_window}", file=sys.stderr)

                vcf_dict = get_softclip_indel_vcf(read, refFasta, row['Pos'], args.mutation_search_window, args.min_softclip_length)

                if (vcf_dict is None or 
                    (vcf_dict['distance'] > args.max_mutation_distance and vcf_dict['distance2'] > args.max_mutation_distance)): 

                    # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                    if unevaluable_read_log:
                        print(f"{str(read)}\t{vcf_dict}", file=unevaluable_read_log)

                    continue

            # read spans start and end and has no softclips    
            elif read.reference_start < row['End'] and read.reference_end > row['Start'] and cigar[0][0] == 0 and cigar[-1][0] == 0:

                if args.verbose:
                    print(f"\tThis read is reference {read.query_name}, {read.cigarstring}, from {row['Chromosome']}:{row['Start']-args.target_window}-{row['End']+args.target_window}", file=sys.stderr)

                vcf_dict = {
                    'read': read.query_name,
                    'chrom': row['Chromosome'],
                    'pos': row['Start'],
                    'distance': '.',
                    'chrom2': row['Chromosome'],
                    'pos2': row['End'],
                    'distance2': '.',
                    'ref': '.',
                    'alt': '.',
                    'strands': '.',
                    'info': '.',
                    'type': 'REF'
                }

            # If read doesnt meet any of the criteria then skip it.
            else:
                # print abbreviated read info (name, cigar, mapping info, sequence) to unevaluable read log
                # if unevaluable_read_log:
                #     print(f"{str(read)}", file=unevaluable_read_log)
                continue

            vcf_dict['info'] = make_info_string(vcf_dict['info']) if vcf_dict['info']!='.' else '.'
            
            # Add indel info to dataframe
            readaln = pd.concat([readaln, pd.DataFrame([vcf_dict])], ignore_index=True)

        if args.verbose:
            print("\tSorting indels", file=sys.stderr)

        # Group by read and process results
        readaln = readaln.sort_values(by=['read','chrom','pos','distance','chrom2','pos2','distance2','strands','ref','alt','type'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
        indelcounts = readaln.groupby(['chrom','pos','distance','chrom2','pos2','distance2','strands','ref','alt','type'],dropna=False).size().reset_index(name='counts')

        # Recast as int type, allowing for NA values
        indelcounts['pos'] = indelcounts['pos'].astype(pd.Int64Dtype())
        indelcounts['pos2'] = indelcounts['pos2'].astype(pd.Int64Dtype())

        if args.verbose:
            print("\tGetting control counts", file=sys.stderr)

        # Get counts from control BAM
        if len(indelcounts) > 0:
            # Add pam positions to the df
            if row['Pos'] is not pd.NA:
                indelcounts['Positions'] = [row['Pos']] * len(indelcounts)
                indelcounts['Distance'] = indelcounts.apply(
                    lambda r: min(
                        [abs(r['pos'] - int(x)) for x in r['Positions']] + 
                        [abs(r['pos'] + len(r['ref']) - 1 - int(x)) for x in r['Positions']]
                    ) if pd.notna(r['pos']) and r['Positions'] else pd.NA,
                    axis=1
                )

            else:
                indelcounts['Positions'] = pd.NA
                indelcounts['Distance'] = pd.NA

            indelcounts = add_normal_counts(indelcounts, [ x for x in control_bamfile.fetch(contig=row['Chromosome'], start=row['Start']-args.target_window, end=row['End']+args.target_window) ], refFasta, window=args.max_mutation_distance)
        
        else:
            indelcounts['control_alt_counts'] = 0
            indelcounts['control_total_counts'] = 0
            indelcounts['Positions'] = pd.NA
            indelcounts['Distance'] = pd.NA

        # Apply filters
        indelcounts = indelcounts[(indelcounts['counts'] >= args.min_coverage) | (indelcounts['ref']=='.')]
        indelcounts = indelcounts[(indelcounts['control_alt_counts'] <= args.max_in_control) | (indelcounts['ref']=='.')]
        indelcounts = indelcounts[(indelcounts['Distance'] <= args.max_mutation_distance) | (indelcounts['ref']=='.')]
        indelcounts = indelcounts[~indelcounts['alt'].str.contains('N')]

        indels = indelcounts

        # # ORIGINAL INDEL CALCULATION
        # indels = get_indels(bam=expsamfile,controlbam=consamfile,chr=row['Chromosome'],start=row['Start'],end=row['End'],
        #                     fasta=refFasta,window=args.window,distance=args.distance,pam_positions=row['Pos'],minreads=args.minreads,maxcontrol=args.maxcontrol,verbose=args.verbose)
        
        # Process indel results
        total_reads, indel_reads, control_total_reads, control_indel_reads = 0, 0, 0, 0
        indel_fraction, control_indel_fraction = 0, 0
        indel_keys, bnd_keys = '.', '.'
        bnds = []

        if len(indels) > 0:
            total_reads = sum(indels['counts'])
            indel_reads = sum(indels[indels['type']!='REF']['counts'])
            control_indel_reads = int(indels['control_alt_counts'].mean())
            control_total_reads = int(indels['control_total_counts'].mean())

            # Separate BNDs and indels
            bnds = indels[indels['type']=='BND'].copy()
            indels = indels[indels['type']=='INDEL'].copy()

            if len(indels) > 0:
                indels['Key'] = indels.apply(lambda r: f"{r['chrom']}|{r['pos']}|{r['chrom2']}|{r['pos2']}|{r['ref']}|{r['alt']}|{r['distance']}|{r['distance2']}|{r['strands']}|{r['counts']}|{r['control_alt_counts']}", axis=1)

            if len(bnds) > 0:
                bnds['Key'] = bnds.apply(lambda r: f"{r['chrom']}|{r['pos']}|{r['chrom2']}|{r['pos2']}|{r['ref']}|{r['alt']}|{r['distance']}|{r['distance2']}|{r['strands']}|{r['counts']}|{r['control_alt_counts']}", axis=1)

        # Calculate fractions and prepare output
        indel_fraction = round(indel_reads/total_reads,4) if total_reads > 0 else 0
        control_indel_fraction = round(control_indel_reads/control_total_reads,4) if control_total_reads > 0 else 0
        indel_keys = ';'.join(indels['Key'].tolist()) if len(indels) > 0 else '.'
        bnd_keys = ';'.join(bnds['Key'].tolist()) if len(bnds) > 0 else '.'
        ontarget = row['Ontarget']
        positions = row['Pos']        
        offtargetsites = make_info_string(row['Info']) if row['Info'] else '.'

        # CRISPR PREDICTION (if enabled)
        crispr_predicted_reads = 0
        crispr_prediction_fraction = 0.0
        crispr_prediction_probability = 0.0
        if args.enable_crispr_prediction and crispr_model:
            if args.verbose:
                print(f"  Predicting CRISPR reads for {row['Chromosome']}:{row['Start']}-{row['End']} (total_reads: {total_reads})", file=sys.stderr)
            
            crispr_predicted_reads, crispr_prediction_probability = predict_reads_at_position(
                edited_bamfile, row['Chromosome'], row['Start'], row['End'], positions, 
                crispr_model, refFasta, is_on_target=ontarget, control_bam=control_bamfile, threshold=args.crispr_threshold
            )
            crispr_prediction_fraction = round(crispr_predicted_reads/total_reads, 4) if total_reads > 0 else 0.0
            
            if args.verbose:
                print(f"    Predicted {crispr_predicted_reads}/{total_reads} reads as CRISPR-related (≥{args.crispr_threshold})", file=sys.stderr)
                print(f"    Average prediction probability: {crispr_prediction_probability:.1f}%", file=sys.stderr)

        # Filter off-target false positives if enabled
        if args.filter_off_target_fp:
            if ontarget == 0 and indel_reads > 0 and crispr_predicted_reads == 0:
                if args.verbose:
                    print(f"  Filtering FP off-target site {row['Chromosome']}:{row['Start']}-{row['End']} (indel_reads: {indel_reads}, predicted_reads: {crispr_predicted_reads})", file=sys.stderr)
                
                if fp_log:
                    log_fields = [
                        row['Chromosome'], row['Start'], row['End'], 
                        ';'.join([str(x) for x in row['Pos']]), total_reads, indel_reads, indel_fraction,
                        control_total_reads, control_indel_reads, control_indel_fraction,
                        len(indels), indel_keys, len(bnds), bnd_keys, offtargetsites, ontarget
                    ]
                    if args.enable_crispr_prediction:
                        log_fields.extend([crispr_predicted_reads, crispr_prediction_fraction, round(crispr_prediction_probability, 1)])
                    fp_log.write("\t".join([str(field) for field in log_fields]) + "\n")
                continue
        
        # Output results
        output_fields = [
            row['Chromosome'], row['Start'], row['End'], 
            ';'.join([str(x) for x in row['Pos']]), total_reads, indel_reads, indel_fraction,
            control_total_reads, control_indel_reads, control_indel_fraction,
            len(indels), indel_keys, len(bnds), bnd_keys, offtargetsites, ontarget
        ]
        if args.enable_crispr_prediction:
            output_fields.extend([crispr_predicted_reads, crispr_prediction_fraction, round(crispr_prediction_probability, 1)])

        print("\t".join([str(field) for field in output_fields]), file=fp, flush=True)

    # ========================================================================
    # STEP 6: Cleanup
    # ========================================================================
    
    fp.close()

    if fp_log:
        fp_log.close()

    if unevaluable_read_log:
        unevaluable_read_log.close()

    edited_bamfile.close()
    control_bamfile.close()
    refFasta.close()

if __name__ == "__main__":
    main()
