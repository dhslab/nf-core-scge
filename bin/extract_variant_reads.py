#!/usr/local/bin/python3

from __future__ import division
import pysam
import biotite.sequence.align as align
import biotite.sequence as seq
import argparse
import re
import string
import csv, sys
import scipy.stats as stats
import pandas as pd, pyranges as pr

# reverse complement function
def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

# make cigar tuples from a cigar string
def make_cigar_tuples(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_tuples = [(cigar_dict[operation],int(length)) for length, operation in cigar_tuples]
    return cigar_tuples

def cigar_summary(cigar):
    cigar_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    # use regex to parse cigar string into tuples
    cigar_tuples = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    cigar_sum = { 'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0 }
    for length, operation in cigar_tuples:
        cigar_sum[operation] += int(length)

    return cigar_sum

def make_aligned_pairs(sup_align,refseq):
    """
    Converts a supplementary alignment string into a list of aligned pairs.

    Args:
    supplementary_alignment (str): Supplementary alignment string (e.g., 'chr12,107092390,-,110M41S,60,0').

    Returns:
    List of tuples: Aligned pairs (e.g., [(0, 168391768, 'C'), (1, 168391769, 'C'), ...]).
    """
    # Split the supplementary alignment string into components
    ref_name, ref_start, strand, cigar, _, _ = supplementary_alignment.split(',')
    
    # Convert positions and CIGAR to appropriate types
    ref_start = int(ref_start)
    cigar_ops = parse_cigar(cigar)
    
    # Initialize variables
    aligned_pairs = []
    read_pos = 0
    ref_pos = ref_start
        
    for op, length in cigar_ops:
        if op == 'M':  # Match/mismatch
            for i in range(length):
                aligned_pairs.append((read_pos, ref_pos, refseq[ref_pos - ref_start]))
                read_pos += 1
                ref_pos += 1
        elif op == 'I':  # Insertion
            for i in range(length):
                aligned_pairs.append((read_pos, None, refseq[ref_pos - ref_start]))
                read_pos += 1
        elif op == 'D':  # Deletion
            for i in range(length):
                aligned_pairs.append((None, ref_pos, refseq[ref_pos - ref_start]))
                ref_pos += 1
        elif op == 'S':  # Soft clipping
            read_pos += length
        # Skipping 'H' (hard clipping) and 'N' (skipped region) as they do not appear in the read

    return aligned_pairs


def indels_from_aligned_pairs(pairs,readseq):

    # remove leading soft clips
    while pairs[0][0] is None:
        pairs = pairs[1:]

    # remove training soft clips
    while pairs[-1][0] is None:
        pairs = pairs[:-1]

    pairs = [ (x[0],x[1],x[2].upper(),None) for x in pairs ]
    for i in range(len(pairs)):
        if pairs[i][0] is not None:
            pairs[i] = (pairs[i][0],readseq[pairs[i][0]],pairs[i][1],pairs[i][2])

    variant_start_index = j-1
    variant_end_index = j-1
    j = 0
    while j < len(pairs):
        if pairs[j][0] is not None and pairs[j][2] is not None:
            j+=1
            continue

        if pairs[j][0] is None or pairs[j][2] is None: # indel
            variant_start_index = j-1
            variant_end_index = len(pairs)
            break

        j+=1

    # construct indel tuple
    indel = (pairs[variant_start_index][2],''.join([ x[3] if x[3] else '' for x in pairs[variant_start_index:variant_end_index]]),''.join([ x[1] if x[1] else '' for x in pairs[variant_start_index:variant_end_index]]))
    
    return indel

# function to count the number of reads that support an indel or BND event
def add_normal_counts(df, reads, fasta, maxreads=0,handicap=5,window=25,flank=150):

    # make ref and alt sequences for each indel/BND
    df['refseq'] = ''
    df['altseq'] = ''
    for it, row in df.iterrows():
        df.at[it,'refseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']+len(row['ref'])+flank)
        if row['type'] in ['INDEL','DEL','INS']:
            df.at[it,'altseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + row['alt'] + fasta.fetch(row['chrom'],row['pos']+len(row['ref'])-1,row['pos']+len(row['ref'])+flank)
        elif row['type'] == 'BND':
            df.at[it,'altseq'] = fasta.fetch(row['chrom'],row['pos']-1-flank,row['pos']-1) + (fasta.fetch(row['chrom2'],row['pos2']-1,row['pos2']+flank) if row['strands']=='++' else revcomp(fasta.fetch(row['chrom2'],row['pos2']-1-flank,row['pos2'])))

    df['control_alt_counts'] = 0
    df['control_total_counts'] = 0

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
        if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == len(read.query_sequence):
            continue
        
        for it, row in df.iterrows():

            if read.cigartuples[0][0]==0 and read.cigartuples[-1][0]==0:
                if (len(row['ref']) > len(row['alt']) and 'D' not in read.cigarstring) or (len(row['alt']) > len(row['ref']) and 'I' not in read.cigarstring):
                    continue

            if row['type'] == 'BND' and not read.has_tag('SA'):
                continue

            refseq = row['refseq']
            altseq = row['altseq']

            ref_align = align.align_optimal(seq.NucleotideSequence(row['refseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]
        
            alt_align = align.align_optimal(seq.NucleotideSequence(row['altseq']),seq.NucleotideSequence(read.query_sequence),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0]

            if alt_align.score > 0 and alt_align.score - handicap > ref_align.score:
                ref_align_cigar_sum = cigar_summary(align.write_alignment_to_cigar(ref_align))
                alt_align_cigar_sum = cigar_summary(align.write_alignment_to_cigar(alt_align))

                if ((len(row['ref'])-len(row['alt']) > 0 and ref_align_cigar_sum['D']==len(row['ref'])-len(row['alt'])) or 
                (len(row['alt'])-len(row['ref']) > 0 and ref_align_cigar_sum['I']==len(row['alt'])-len(row['ref']))):
                    df.at[it,'control_alt_counts'] += 1
                    # print to stderr: found control alt count
                    print(f"\tFound control alt count for {row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}:{read.query_name}", file=sys.stderr)

    df['control_total_counts'] = len(total_reads)

    return df.copy()


# function to get indels for one amplicon from bam file
def get_indels(bam,controlbam,chr,start,end,fasta,window=100,distance=25,pam_positions=pd.NA,minSecMapQual=10,maxNM=5,svDistanceThreshold=100000,saAlignmentTolerance=5,minSoftClipLength=5,deletionDistanceThreshold=1000,minreads=1,maxcontrol=0,strict=False,verbose=False):

    # window == flanking sequence to add to read search
    # distance == max distance from start/end to consider an indel, BND, or call a read wild-type.

    cigarVarDict = {'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}

    readaln = pd.DataFrame(columns=['read','chrom','pos','chrom2','pos2','ref','alt','strands','type'])
    
    # if verbose, print to stderr: getting reads that align within window of start and end
    if verbose:
        print(f"\tGetting reads that align within window of {chr}:{start}-{end}", file=sys.stderr)

    regionStart = max(start - svDistanceThreshold,1)
    regionEnd = min(end + svDistanceThreshold,fasta.get_reference_length(chr))
    regionSeq = fasta.fetch(chr,regionStart-1,regionEnd)

    # get reads that align within window of start and end
    for read in bam.fetch(chr, start-window, end+window, multiple_iterators = True):

        # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
        if read.is_mapped is False or \
            read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.mapping_quality == 0:
            continue

        cigar = read.cigartuples # get cigar info
        mate_cigar = make_cigar_tuples(read.get_tag('MC')) if read.has_tag('MC') else None

        # Quality filter: Get mismatches from NM tag and subtract deleted/skipped bases from cigar
        # and skip if too many mismatches
        if read.has_tag('NM') and int(read.get_tag('NM')) > maxNM:
            if (int(read.get_tag('NM')) - sum([x[1] for x in cigar if x[0] == 2 or x[0]==3])) > maxNM:
                continue

        read_strand = "+"
        if read.is_reverse is True:
            read_strand = "-"

        leftSoftClip = 0
        rightSoftClip = 0
        read_reference_start = read.reference_start + 1
        read_reference_end = read.reference_end
        read_next_reference_start = read.next_reference_start + 1
        read_next_reference_end = None if mate_cigar is None else read_next_reference_start + sum([x[1] for x in mate_cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])
        
        read_query_start = 0
        read_query_end = len(read.query_sequence)

        # setup indelinfo dictionary
        indelinfo = {'read':read.query_name,'chrom':'','pos':None,'chrom2':None,'pos2':None,'strands': '','ref':'','alt':'','type':''}

        # Get left soft clip length
        if cigar[0][0] == 4:
            leftSoftClip = cigar[0][1]

        # Get right soft clip length
        if cigar[-1][0] == 4:
            rightSoftClip = cigar[-1][1]

        # 
        # This adjusts alignment cigars for soft clips and supplementary alignments. 
        # It extends the cigar if the softclips map nearby and therefore the read has a deletion.
        #

        # Process left soft clip, if present
        if leftSoftClip > rightSoftClip:

            # first check for local supplementary alignments 
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None and int(sMq)>=minSecMapQual and
                int(sNm)<maxNM and sChr == read.reference_name and
                sStrand == read_strand and int(sPos) < read_reference_start and
                abs(int(sPos) - read_reference_start) < svDistanceThreshold):
                
                sPos = int(sPos)
                sCigarTuples = make_cigar_tuples(sCigar)

                if sCigarTuples[-1][0] >= 4 and sCigarTuples[0][0] == 0 and sCigarTuples[0][1] == cigar[0][1]:
                    cigar = sCigarTuples[:-1] + [(2,read_reference_start-sPos-1)] + cigar[1:]

                else:
                    # if there is a supplementary alignment, then realign the entire read to the reference sequence
                    cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[sPos-regionStart:read_reference_end-regionStart+1]), #seq.NucleotideSequence(fasta.fetch(chr,sPos-1,read_reference_end)),
                                                                                seq.NucleotideSequence(read.query_sequence),
                                                                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                                                gap_penalty=(-10,-1),local=False,max_number=1)[0]))
                read_reference_start = sPos

            else:
                # if the left soft clip is long enough, see if it supports a deletion within window
                # this performs an exact match and so isnt ideal, but should be good enough for now
                if leftSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[:leftSoftClip]
                    refSeq = regionSeq[read_reference_start - regionStart - window + 1:read_reference_start - regionStart +  1] #fasta.fetch(chr,read_reference_start - window,read_reference_start)
                    # find the position of the left-most occurance of clippedSeq in refSeq
                    leftClipPos = refSeq.rfind(clippedSeq)
                    if leftClipPos != -1:
                        delSeq = refSeq[leftClipPos + len(clippedSeq)-1:]
                        cigar = [(0,len(clippedSeq))] + [(2,len(delSeq))] + cigar[1:]
                        read_reference_start = read_reference_end - sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3]) + 1

        # Process right soft clip, if present
        if rightSoftClip > leftSoftClip:

            # first check for local supplementary alignments 
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None):
                sPos = int(sPos)
                sCigarTuples = make_cigar_tuples(sCigar)

                if (int(sMq)>=minSecMapQual and sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and
                    int(sNm)<maxNM and sChr == read.reference_name and 
                    sStrand == read_strand and int(sPos) > read_reference_end and 
                    int(sPos) - read_reference_end < svDistanceThreshold):
                
                    if sCigarTuples[0][0] >= 4 and sCigarTuples[-1][0] == 0 and sCigarTuples[-1][1] == cigar[-1][1]:
                        cigar = cigar[:-1] + [(2,read_reference_end-sPos)] + sCigarTuples[1:]

                    else:
                        cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(regionSeq[read_reference_start-regionStart:sPos-regionStart+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3])]), #fasta.fetch(chr,read_reference_start-1,sPos-1+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]))),
                                        seq.NucleotideSequence(read.query_sequence),
                                        matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                        gap_penalty=(-10,-1),local=False,max_number=1)[0]))

                read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])

            else:
                # if the left soft clip is long enough, see if it supports a deletion within deletionDistanceThreshold
                # this performs an exact match and so isnt ideal, but should be good enough for now
                if rightSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[-rightSoftClip:]
                    refSeq = regionSeq[read_reference_end-regionStart:read_reference_end-regionStart + window] #fasta.fetch(chr,read_reference_end,read_reference_end + window)
                    # find the position of the left-most occurance of clippedSeq in refSeq
                    rightClipPos = refSeq.find(clippedSeq)
                    if rightClipPos != -1:
                        delSeq = refSeq[0:rightClipPos]
                        cigar = cigar[:-1] + [(2,len(delSeq))] + [(0,len(clippedSeq))]
                        read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])
        
        # Adjusted read has to be within distance of the start/end positions to consider
        if read_reference_start - leftSoftClip >= end + distance or read_reference_end + rightSoftClip <= start - distance:
            continue

        # remove leading or trailing soft clips
        if cigar[0][0] >= 4:
            read_query_start += cigar[0][1]
            del cigar[0]
        
        if cigar[-1][0] >= 4:
            read_query_end -= cigar[-1][1]
            del cigar[-1]

        # there are multiple cigar operations bookened by matches, the process an indel                
        if len(cigar) > 1:

            # adjust alignment start for first aligned block and remove it
            if cigar[0][0] == 0:
                read_reference_start += cigar[0][1]# sets to position before event
                read_query_start = read_query_start + cigar[0][1] 
                del cigar[0]

            # adjust alignment end for last aligned block and remove it
            if cigar[-1][0] == 0:
                read_reference_end = read_reference_end - cigar[-1][1]
                read_query_end = read_query_end - cigar[-1][1]
                del cigar[-1]

            # If there's only one cigar operation left, it's an insertion or deletion
            if len(cigar) == 1:

                # theres one event, so record the coordinates
                indelinfo['chrom'] = read.reference_name
                indelinfo['pos'] = read_reference_start
                indelinfo['strands'] = '++'

                # if insertion
                if cigar[0][0] == 1:           
                    indelinfo['ref'] = regionSeq[read_reference_start-regionStart-1:read_reference_start-regionStart] #fasta.fetch(chr,read_reference_start-1-1,read_reference_start-1) # 0 based and position before the isertion
                    indelinfo['alt'] = read.query_sequence[read_query_start-1:read_query_start+cigar[0][1]]
                    indelinfo['type'] = 'INDEL'

                # if deletion
                elif cigar[0][0] == 2 or cigar[0][0] == 3:
                    indelinfo['ref'] = regionSeq[read_reference_start-regionStart:read_reference_start-regionStart + cigar[0][1] + 1] #fasta.fetch(chr,read_reference_start-1,read_reference_start+cigar[0][1])
                    indelinfo['alt'] = indelinfo['ref'][0]
                    indelinfo['type'] = 'INDEL'

            # if its a complex indel
            else:
                indelinfo['chrom'] = read.reference_name
                indelinfo['pos'] = read_reference_start
                indelinfo['strands'] = '++'
                indelinfo['type'] = 'INDEL'

                # iterate through cigar and get indels
                varlen = sum([ x[1] for x in cigar ])
                indelinfo['ref'] = regionSeq[read_reference_start-regionStart:read_reference_start-regionStart+sum([ x[1] for x in cigar ])+1] #fasta.fetch(chr,read_reference_start-1,read_reference_start+sum([ x[1] for x in cigar ]))
                indelinfo['alt'] = read.query_sequence[read_query_start+1:read_query_start+sum([x[1] if x[0] == 1 or x[0] == 0 else 0 for x in cigar])]

        # if the read is chimeric and has a supplementary alignment near the start/end, process the chimeric read
        elif read.has_tag('SA') is True and (leftSoftClip > 0 or rightSoftClip > 0):
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            sCigarTuples = make_cigar_tuples(sCigar)
            sPos = int(sPos)
            sEnd = sPos + sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]) - 1
            sLeftSoftClip = sCigarTuples[0][1] if sCigarTuples[0][0] >= 4 else 0
            sRightSoftClip = sCigarTuples[-1][1] if sCigarTuples[-1][0] >= 4 else 0

            pSeq = read.query_alignment_sequence
            sSeq = read.query_sequence[sLeftSoftClip:sLeftSoftClip+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 1])] if sStrand == read_strand else revcomp(read.query_sequence[sLeftSoftClip:sLeftSoftClip+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 1])])

            # correct for overlapping alignments
            if len(pSeq)+len(sSeq) > len(read.query_sequence):
                trimLen = len(pSeq)+len(sSeq) - len(read.query_sequence)
                if sLeftSoftClip > sRightSoftClip:
                    sPos += trimLen
                    sSeq = sSeq[trimLen:]

                elif sRightSoftClip > sLeftSoftClip:
                    sEnd -= trimLen
                    sSeq = sSeq[:-trimLen]

            if (int(sMq) >= minSecMapQual and int(sNm) <= maxNM and int(sMq)>0 and
                (sChr != read.reference_name or abs(int(sPos) - read_reference_start) >= svDistanceThreshold) or sStrand != read_strand or 
                (rightSoftClip > 0 and sEnd < read_reference_start) or (leftSoftClip > 0 and sPos > read_reference_end)):

                indelinfo['chrom'] = read.reference_name
                indelinfo['chrom2'] = sChr
                indelinfo['strands'] = '++'if sStrand == read_strand else '+-'
                indelinfo['type'] = 'BND'

                if indelinfo['strands'] == '++':
                    # if --->...> ...>---> or <---<... <...<---
                    if rightSoftClip > leftSoftClip and sLeftSoftClip > sRightSoftClip:
                        indelinfo['pos'] = read_reference_end
                        indelinfo['pos2'] = sPos

                    # if ...>---> --->...> or <...<--- <---<...
                    elif leftSoftClip > rightSoftClip and sRightSoftClip > sLeftSoftClip:
                        indelinfo['pos'] = read_reference_start - 1 
                        indelinfo['pos2'] = sEnd + 1

                    else:
                        print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                        if strict:
                            exit(1)
                        else:
                            continue

                elif indelinfo['strands'] == '+-':
                    # if --->...> <---<... or <---<... --->...>
                    if rightSoftClip > leftSoftClip and sRightSoftClip > sLeftSoftClip:
                        indelinfo['pos'] = read_reference_end
                        indelinfo['pos2'] = sEnd

                    # if ...>---> <...<--- or <...<--- ...>--->
                    elif leftSoftClip > rightSoftClip and sLeftSoftClip > sRightSoftClip:
                        indelinfo['pos'] = read_reference_start - 1
                        indelinfo['pos2'] = sPos

                    else:
                        print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                        if strict:
                            exit(1)
                        else:
                            continue
                
                else:
                    print(f"Error: no proper soft clip orientation found: {read.query_name}: {str(read)}", file=sys.stderr)
                    if strict:
                        exit(1)
                    else:
                        continue
            
            elif read_reference_start < end and read_reference_end > start:
                indelinfo['chrom'] = chr
                indelinfo['pos'] = start
                indelinfo['ref'] = '.'
                indelinfo['alt'] = '.'
                indelinfo['type'] = 'REF'

            else:
                continue
                    
        elif read_reference_start < end and read_reference_end > start:
            indelinfo['chrom'] = chr
            indelinfo['pos'] = start
            indelinfo['ref'] = '.'
            indelinfo['alt'] = '.'
            indelinfo['type'] = 'REF'

        else:
           continue

        # add indel info to dataframe
        readaln = pd.concat([readaln, pd.DataFrame([indelinfo])], ignore_index=True)

    # if verbose, sorting indels
    if verbose:
        print("\tSorting indels", file=sys.stderr)

    # group by read and sort by cigar, then sv and get the first value in each group
    readaln = readaln.sort_values(by=['read','chrom','pos','chrom2','pos2','strands','ref','alt','type'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
    indelcounts = readaln.groupby(['chrom','pos','chrom2','pos2','strands','ref','alt','type'],dropna=False).size().reset_index(name='counts')

    # need to recast as int type, but allow for NA values
    indelcounts['pos'] = indelcounts['pos'].astype(pd.Int64Dtype())
    indelcounts['pos2'] = indelcounts['pos2'].astype(pd.Int64Dtype())

    # if verbose, get control counts
    if verbose:
        print("\tGetting control counts", file=sys.stderr)

    # get counts of reads in the controlbam for each indel/BND
    if len(indelcounts) > 0:
        # add pam positions to the df
        if pam_positions is not pd.NA:
            indelcounts['Positions'] = [pam_positions] * len(indelcounts)
            indelcounts['Distance'] = indelcounts.apply(lambda r: min([abs(r['pos'] - x) for x in r['Positions']] + [abs(r['pos'] + len(r['ref']) - 1 - x) for x in r['Positions']] ) if r['pos'] is not pd.NA else pd.NA,axis=1)
        else:
            indelcounts['Positions'] = pd.NA
            indelcounts['Distance'] = pd.NA

        indelcounts = add_normal_counts(indelcounts, [ x for x in controlbam.fetch(contig=chr,start=start-window,end=end) ], fasta, window=distance)

    else:
        indelcounts['control_alt_counts'] = 0
        indelcounts['control_total_counts'] = 0
        indelcounts['Positions'] = pd.NA
        indelcounts['Distance'] = pd.NA


    # apply filters
    indelcounts = indelcounts[(indelcounts['counts'] >= minreads) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['control_alt_counts'] <= maxcontrol) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['Distance'] <= distance) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[~indelcounts['alt'].str.contains('N')]

    return(indelcounts)


def main():

    parser = argparse.ArgumentParser(description='Find indels in a bam file at BED coordinates')
    parser.add_argument('-f','--fasta',type=str,default="/storage1/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh_v4.3.6/hg38_PLVM_CD19_CARv4_cd34.fa",help='Reference fasta file')
    parser.add_argument('-w','--window',type=int,default=100,help='Distance between off-target sites for merging intervals.')
    parser.add_argument('-d','--distance',type=int,default=25,help='Window size around off-target site to identify mutations')
    parser.add_argument('-m','--minreads',type=int,default=1,help='Minimum supporting reads to report an indel/bnd event.')
    parser.add_argument('-c','--chromosome',type=str,default=None,help='Chromosome to process.')
    parser.add_argument('-x','--maxcontrol',type=int,default=0,help='Maximum supporting reads to to report an indel/bnd event.')
    parser.add_argument('-s','--strict',action='store_true',help='Exit if a read cannot be properly parsed.')
    # add verbosity argument
    parser.add_argument('-v','--verbose',action='store_true',help='Print verbose output')

    # add outfile argument -o or --outfile
    parser.add_argument('-o','--outfile',type=str,help='Output file (optional)')
    parser.add_argument('bed',help='Off-target coordinate file')
    parser.add_argument('expbamfile',help='BAM file')
    parser.add_argument('conbamfile',help='BAM file')

    args = parser.parse_args()

    # TODO: 
    # Nihdi: We need to update this code to take the file and return 2 objects:
    # 1. Take the positions and make intervals that are chr,(pos-distup),(pos+distdown) and then merge them into new list of non-redunant intervals (might use pyranges for this)
    # 2. The individual entries need to be retained, though, along with the information about the off-target (the exact position, the mismatches, etc)

    # verbose to stderr: processing input file
    if args.verbose:
        print("Processing input file", file=sys.stderr)

    # open off-target file and create ranges to search
    bedDf = pd.read_csv(args.bed)
    bedDf['Pos'] = bedDf['Start']
    bedDf['End'] = bedDf['Start']
    bedDf['Start'] = bedDf['Start'] - 1
    bedDf['Info'] = bedDf.apply(lambda row: f"{row['Source']},{row['DNA Sequence']},{row['PAM']},{row['Chromosome']},{row['Pos']},{row['Strand Direction']},{row['Mismatch']},{row['Bulge Type']},{row['Bulge Size']}", axis=1)
    bedDf['Ontarget'] = bedDf['On_target']

    # make pyranges object
    bedPr = pr.PyRanges(bedDf[['Chromosome','Start','End','Pos','Info','Ontarget']])

    # cluster the intervals
    bedPr = bedPr.cluster(slack=args.window)
    mergedBedPr = bedPr.merge(by='Cluster',strand=False,slack=args.window)
    mergedBedDf = mergedBedPr.df.join(bedPr.df.groupby('Cluster')['Pos'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')
    mergedBedDf = mergedBedDf.join(bedPr.df.groupby('Cluster')['Info'].agg(list).reset_index().set_index('Cluster'),on='Cluster',how='left')
    mergedBedDf = mergedBedDf.join(bedPr.df.groupby('Cluster')['Ontarget'].agg('first').reset_index().set_index('Cluster'),on='Cluster',how='left')

    if args.chromosome is not None:
        print("Processing chromosome", args.chromosome, file=sys.stderr)
        mergedBedDf = mergedBedDf[mergedBedDf['Chromosome'] == args.chromosome]

    # verbose: say done
    if args.verbose:
        print("Done processing input file", file=sys.stderr)

    # open bam file(s)
    expsamfile = pysam.AlignmentFile(args.expbamfile,"rc",reference_filename=args.fasta)
    consamfile = pysam.AlignmentFile(args.conbamfile,"rc",reference_filename=args.fasta)
    # open fasta file
    refFasta = pysam.FastaFile(args.fasta)

    # print to outfile or stdout
    if args.outfile:
        sys.stdout = open(args.outfile, 'w')

    print("\t".join('chrom start end pam_positions total_reads indel_reads indel_fraction control_reads control_indel_reads control_indel_fraction indel_count indel_info bnd_count bnd_info target_info is_target'.split(' ')),flush=True)

    # iterate over mergedBedPr intervals:
    for index, row in mergedBedDf.iterrows():

        # verbose to stderr: processing interval
        if args.verbose:
            print(f"Processing interval {row['Chromosome']}:{row['Start']}-{row['End']}", file=sys.stderr)

        indels = get_indels(bam=expsamfile,controlbam=consamfile,chr=row['Chromosome'],start=row['Start'],end=row['End'],
                            fasta=refFasta,window=args.window,distance=args.distance,pam_positions=row['Pos'],minreads=args.minreads,maxcontrol=args.maxcontrol,verbose=args.verbose)
        
        total_reads, indel_reads, control_total_reads, control_indel_reads = 0, 0, 0, 0
        indel_fraction, control_indel_fraction = 0, 0
        indel_keys, bnd_keys = '.', '.'
        bnds = []

        if len(indels) > 0:
            # get total reads in this region
            total_reads = sum(indels['counts'])
            indel_reads = sum(indels[indels['type']!='REF']['counts'])
            control_indel_reads = int(indels['control_alt_counts'].mean())
            control_total_reads = int(indels['control_total_counts'].mean())

            # separate BNDs and indels
            bnds = indels[indels['type']=='BND'].copy()
            indels = indels[indels['type']=='INDEL'].copy()

            if len(indels) > 0:
                indels['Key'] = indels.apply(lambda r: f"{r['chrom']}:{r['pos']}:{r['ref']}:{r['alt']}:{r['counts']}:{r['control_alt_counts']}:{r['Distance']}", axis=1)

            if len(bnds) > 0:
                bnds['Key'] = bnds.apply(lambda r: f"{r['chrom']}:{r['pos']}:{r['chrom2']}:{r['pos2']}:{r['strands']}:{r['counts']}:{r['control_alt_counts']}:{r['Distance']}", axis=1)
    
        # print results to tab-delimited file
        indel_fraction = round(indel_reads/total_reads,4) if total_reads > 0 else 0
        control_indel_fraction = round(control_indel_reads/control_total_reads,4) if control_total_reads > 0 else 0
        indel_keys = ';'.join(indels['Key'].tolist()) if len(indels) > 0 else '.'
        bnd_keys = ';'.join(bnds['Key'].tolist()) if len(bnds) > 0 else '.'
        ontarget = row['Ontarget']
        offtargetsites = ';'.join(row['Info']) if len(row['Info']) > 0 else '.'

        print(f"{row['Chromosome']}\t{row['Start']}\t{row['End']}\t{';'.join([ str(x) for x in row['Pos'] ])}\t{total_reads}\t{indel_reads}\t{indel_fraction}\t{control_total_reads}\t{control_indel_reads}\t{control_indel_fraction}\t{len(indels)}\t{indel_keys}\t{len(bnds)}\t{bnd_keys}\t{offtargetsites}\t{ontarget}",flush=True)


    # close stdout or outfile
    if args.outfile:
        sys.stdout.close()

    expsamfile.close()
    consamfile.close()
    refFasta.close()


if __name__ == "__main__":
    main()