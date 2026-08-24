#!/usr/bin/env python3
"""
wgs_tag_reads.py — per-read XC tags for the WGS arm, for live IGV review.

The ECS arm has had this since the read-level-tags work: every read in the target windows
carries an `XC` tag naming how the caller classified it, so a reviewer can colour the pileup
(Color alignments by -> tag -> XC) and see *why* a site was called instead of taking the
caller's word for it. The WGS arm had only rendered PNGs, which cannot be interrogated.

This closes that gap, and deliberately reuses the ECS machinery — `merge_windows`,
`write_tagged_bam`, `record_read_tag`, `read_tag_key` are imported from find_edited_reads.py
rather than reimplemented — so **one IGV colour scheme works for both arms**. Only the per-read
classification is WGS-specific, because the WGS caller works from pileup shape at a candidate
locus rather than from per-read VCF records.

Why not reuse features.py:read_records()
---------------------------------------
It returns `{'spans','indel','mapq','softclip'}` with the read identity discarded, and it
`continue`s past duplicate/secondary/supplementary reads — so the `Skipped_*` classes, which are
the most useful thing a reviewer can see, never reach it. It is also the version-pinned model
input path. A separate read walk is both easier and safer than perturbing that.

Tag vocabulary (identical to the ECS arm)
-----------------------------------------
  Edited_Deletion_<N>bp / Edited_Insertion_<N>bp   read carries an indel at the cut
  Edited_SoftClip                                  soft-clip edge at the cut (realignment)
  Unedited_WT                                      spans the cut, no indel: the reference class
  Skipped_Duplicate / Skipped_LowMapQ /
  Skipped_Mismatches / Skipped_NoSpan              excluded, and visibly so
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from find_edited_reads import (read_tag_key, record_read_tag,  # noqa: E402
                               merge_windows, write_tagged_bam,
                               count_mismatches_fast)

# Match the WGS caller's own read filters (features.py read_records / score.py defaults) so the
# tags explain the SAME denominator the features were computed from. If these drift from the
# caller, the BAM stops being an explanation of the call.
MIN_MAPQ = 20
MAX_NM = 4
WINDOW = 150          # bp each side of the cut to emit; matches the ECS target window
CUT_SLACK = 25        # an indel this close to the cut counts as "at" the cut


def classify_wgs_read(read, cut_pos, cut_slack=CUT_SLACK,
                      min_mapq=MIN_MAPQ, max_nm=MAX_NM):
    """One read -> one XC tag string, in the caller's own precedence order."""
    if read.is_duplicate:
        return "Skipped_Duplicate"
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return "Skipped_Unevaluable"
    if (read.mapping_quality or 0) < min_mapq:
        return "Skipped_LowMapQ"
    # NM counts indel bases as well as substitutions, so a genuine 5 bp deletion carries NM>=5
    # and a raw `NM > max_nm` test rejects precisely the reads we are looking for. The ECS
    # caller solves this with count_mismatches_fast(), which subtracts the CIGAR indel lengths
    # to leave true single-base mismatches; reuse it rather than re-deriving it.
    if count_mismatches_fast(read) > max_nm:
        return "Skipped_Mismatches"
    end = read.reference_end or read.reference_start
    if not (read.reference_start <= cut_pos <= end):
        return "Skipped_NoSpan"

    # walk the CIGAR for the largest indel near the cut, plus soft-clip edges
    refpos = read.reference_start
    best = None                     # (op, length)
    softclip_at_cut = False
    for op, ln in (read.cigartuples or []):
        if op in (0, 7, 8):                       # M/=/X consume ref and query
            refpos += ln
        elif op == 2:                             # D consumes ref only
            if abs(refpos - cut_pos) <= cut_slack and (best is None or ln > best[1]):
                best = (2, ln)
            refpos += ln
        elif op == 1:                             # I consumes query only
            if abs(refpos - cut_pos) <= cut_slack and (best is None or ln > best[1]):
                best = (1, ln)
        elif op == 4:                             # S soft-clip; record the ref edge
            if abs(refpos - cut_pos) <= cut_slack:
                softclip_at_cut = True
    if best is not None:
        return ("Edited_Deletion_%dbp" if best[0] == 2 else "Edited_Insertion_%dbp") % best[1]
    if softclip_at_cut:
        return "Edited_SoftClip"
    return "Unedited_WT"


def tag_wgs_bam(bam, sites, out_path, tag_name="XC", window=WINDOW,
                cut_slack=CUT_SLACK, verbose=True):
    """Tag every read in the windows around `sites` and write one indexed BAM.

    `bam`   an open pysam.AlignmentFile (reused from the caller's cache — reopening a CRAM
            per site is the slow way to do this).
    `sites` iterable of (chrom, cut_pos) 0/1-based-agnostic reference positions.
    Returns the number of records written.
    """
    read_tags, windows = {}, []
    for chrom, cut in sites:
        chrom, cut = str(chrom), int(cut)
        lo, hi = max(0, cut - window), cut + window
        windows.append((chrom, lo, hi))
        try:
            fetched = bam.fetch(chrom, lo, hi)
        except (ValueError, KeyError):
            # contig absent from this CRAM (alt contigs are a real case here)
            continue
        for read in fetched:
            # record_read_tag keeps the highest-precedence tag when windows overlap, so a read
            # visited by two nearby candidate sites is still written exactly once.
            record_read_tag(read_tags, read, classify_wgs_read(read, cut, cut_slack))
    if not read_tags:
        if verbose:
            print("  no reads in any candidate window; not writing a WGS tagged BAM",
                  file=sys.stderr)
        return 0
    return write_tagged_bam(bam, out_path, windows, read_tags, tag_name, verbose=verbose)
