#!/usr/bin/python3

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

# function to count the number of reads that support an indel or BND event
def calculate_normal_counts(row, bam, fasta, window=25):

    # subfunctoin to count read hits to a sequence
    def count_read_hits(bam,fasta,chrom,pos,chrom2=None,pos2=None,ref=None,alt=None,orientation=None,handicap=5,window=25,flank=150):

        # note, the handicap is used to adjust the score of the alignment to account for the difference in length between the ref and alt sequences
        # if the normal reads really support the indel, it should have a score that reflects the difference in length between the ref and alt sequences
        refseq = ''
        altseq = ''

        # make ref and alt sequences, first if an indel, next if a BND.
        # note that if its a BND then the refseq is a concatenation of the two sequences from either end
        if ref and alt:
            refseq = fasta.fetch(chrom,pos-1-flank,pos+len(ref)+flank)
            altseq = fasta.fetch(chrom,pos-1-flank,pos-1) + alt + fasta.fetch(chrom,pos+len(ref)-1,pos+len(ref)+flank)
            #handicap = abs(len(ref) - len(alt)) if handicap == 0 else handicap

        elif not chrom2 is pd.NA and not pos2 is pd.NA and not orientation is None:
            handicap = 50 if handicap == 0 else handicap
            refseq = fasta.fetch(chrom,pos-1-flank,pos+flank)
            altseq = fasta.fetch(chrom,pos-1-flank,pos-1)
            if orientation == '++':
                altseq += fasta.fetch(chrom2,pos2-1,pos2+flank)
            elif orientation == '+-':
                altseq += revcomp(fasta.fetch(chrom2,pos2-1-flank,pos2))

        # make ref_count and alt_count counters
        ref_count = set()
        alt_count = set()

        # iterate through reads in bam file for the first position
        for read in bam.fetch(chrom,pos-1-window,pos+window,multiple_iterators=True):
            if not read.is_mapped:
                continue

            # if the read has no indels, add it to the ref_count set
            if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == len(read.query_sequence):
                ref_count.add(read.query_name)
                continue

            ref_align = align.align_optimal(seq.NucleotideSequence(read.query_sequence),seq.NucleotideSequence(refseq),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0].score

            alt_align = align.align_optimal(seq.NucleotideSequence(read.query_sequence),seq.NucleotideSequence(altseq),matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),gap_penalty=(-10,-1),local=True,terminal_penalty=False,max_number=1)[0].score

            if alt_align > 0 and alt_align - handicap > ref_align:
                alt_count.add(read.query_name)
            else:
                ref_count.add(read.query_name)

        return (len(alt_count),len(alt_count)+len(ref_count))

    if row['type'] == 'BND':
        return count_read_hits(
            bam, fasta, window=window,
            chrom=row['chrom'], pos=row['pos'],
            chrom2=row['chrom2'], pos2=row['pos2'],
            orientation=row['strands']
        )
    elif row['type'] == 'INDEL':
        return count_read_hits(
            bam, fasta, window=window,
            chrom=row['chrom'], pos=row['pos'],
            ref=row['ref'], alt=row['alt']
        )
    else:
        #print(f"{row['chrom']}\t{row['pos']}\t{row['type']}")
        control_reads = [ x.query_name for x in bam.fetch(row['chrom'],row['pos']-1-window,row['pos']+window,multiple_iterators=True) ]
        control_reads = set(control_reads)
        return (0,len(control_reads))


# function to get indels for one amplicon from bam file
def get_indels(bam,controlbam,chr,start,end,fasta,window=100,distance=25,pam_positions=pd.NA,minSecMapQual=10,maxNM=5,svDistanceThreshold=100000,saAlignmentTolerance=5,minSoftClipLength=5,deletionDistanceThreshold=1000,minreads=1,maxcontrol=0,strict=False):

    # window == flanking sequence to add to read search
    # distance == max distance from start/end to consider an indel, BND, or call a read wild-type.

    cigarVarDict = {'0':'M','1':'I','2':'D','3':'N','4':'S','5':'H','6':'P','7':'=','8':'X'}

    readaln = pd.DataFrame(columns=['read','chrom','pos','chrom2','pos2','ref','alt','strands','type'])
    
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
        if leftSoftClip > 0:

            # first check for local supplementary alignments 
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None and int(sMq)>=minSecMapQual and 
                int(sNm)<maxNM and sChr == read.reference_name and 
                sStrand == read_strand and int(sPos) < read_reference_start and 
                abs(int(sPos) - read_reference_start) < svDistanceThreshold):

                sPos = int(sPos)

                # if there is a supplementary alignment, then realign the entire read to the reference sequence
                cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(fasta.fetch(chr,sPos-1,read_reference_end)),
                                                                            seq.NucleotideSequence(read.query_sequence),
                                                                            matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                                                            gap_penalty=(-10,-1),local=False,max_number=1)[0]))
                read_reference_start = sPos

            else:
                # if the left soft clip is long enough, see if it supports a deletion within window
                # this performs an exact match and so isnt ideal, but should be good enough for now
                if leftSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[:leftSoftClip]
                    refSeq = fasta.fetch(chr,read_reference_start - window,read_reference_start)
                    # find the position of the left-most occurance of clippedSeq in refSeq
                    leftClipPos = refSeq.rfind(clippedSeq)
                    if leftClipPos != -1:
                        delSeq = refSeq[leftClipPos + len(clippedSeq):]
                        cigar = [(0,len(clippedSeq))] + [(2,len(delSeq))] + cigar[1:]
                        read_reference_start = read_reference_end - sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3]) + 1

        # Process right soft clip, if present
        if rightSoftClip > 0:

            # first check for local supplementary alignments 
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(';')[0].split(',') if read.has_tag('SA') else [None,None,None,None,None,None]

            if (sChr != None and int(sMq)>=minSecMapQual and 
                int(sNm)<maxNM and sChr == read.reference_name and 
                sStrand == read_strand and int(sPos) > read_reference_start and 
                int(sPos) - read_reference_start < svDistanceThreshold):
                
                sPos = int(sPos)
                sCigarTuples = make_cigar_tuples(sCigar)

                cigar = make_cigar_tuples(align.write_alignment_to_cigar(align.align_optimal(seq.NucleotideSequence(fasta.fetch(chr,read_reference_start-1,sPos-1+sum([x[1] for x in sCigarTuples if x[0] == 0 or x[0] == 2 or x[0] == 3]))),seq.NucleotideSequence(read.query_sequence),
                                matrix=align.SubstitutionMatrix.std_nucleotide_matrix(),
                                gap_penalty=(-10,-1),local=False,max_number=1)[0]))

                read_reference_end = read_reference_start + sum([x[1] for x in cigar if x[0] == 0 or x[0] == 2 or x[0] == 3])

            else:
                # if the left soft clip is long enough, see if it supports a deletion within deletionDistanceThreshold
                # this performs an exact match and so isnt ideal, but should be good enough for now
                if rightSoftClip >= minSoftClipLength:
                    clippedSeq = read.query_sequence[-rightSoftClip:]
                    refSeq = fasta.fetch(chr,read_reference_end,read_reference_end + window)
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
                    indelinfo['ref'] = fasta.fetch(chr,read_reference_start-1-1,read_reference_start-1) # 0 based and position before the isertion
                    indelinfo['alt'] = read.query_sequence[read_query_start-1:read_query_start+cigar[0][1]]
                    indelinfo['type'] = 'INDEL'

                # if deletion
                elif cigar[0][0] == 2 or cigar[0][0] == 3:
                    indelinfo['ref'] = fasta.fetch(chr,read_reference_start-1,read_reference_start+cigar[0][1])
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
                indelinfo['ref'] = fasta.fetch(chr,read_reference_start-1,read_reference_start+sum([ x[1] for x in cigar ]))
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
                        print(f"Error: no proper soft clip orientation found: {read.query_name}")
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
                        print(f"Error: no proper soft clip orientation found: {read.query_name}")
                        if strict:
                            exit(1)
                        else:
                            continue
                
                else:
                    print(f"Error: no proper soft clip orientation found: {read.query_name}")
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

    # group by read and sort by cigar, then sv and get the first value in each group
    readaln = readaln.sort_values(by=['read','chrom','pos','chrom2','pos2','strands','ref','alt','type'],key=lambda col: col != '',ascending=False).groupby('read').first().reset_index()
    indelcounts = readaln.groupby(['chrom','pos','chrom2','pos2','strands','ref','alt','type'],dropna=False).size().reset_index(name='counts')

    # need to recast as int type, but allow for NA values
    indelcounts['pos'] = indelcounts['pos'].astype(pd.Int64Dtype())
    indelcounts['pos2'] = indelcounts['pos2'].astype(pd.Int64Dtype())

    # get counts of reads in the controlbam for each indel/BND
    if len(indelcounts) > 0:
        # add pam positions to the df
        if pam_positions is not pd.NA:
            indelcounts['Positions'] = [pam_positions] * len(indelcounts)
            indelcounts['Distance'] = indelcounts.apply(lambda r: min([abs(r['pos'] - x) for x in r['Positions']] + [abs(r['pos'] + len(r['ref']) - 1 - x) for x in r['Positions']] ) if r['pos'] is not pd.NA else pd.NA,axis=1)
        else:
            indelcounts['Positions'] = pd.NA
            indelcounts['Distance'] = pd.NA

        indelcounts[['control_alt_counts','control_total_counts']] = indelcounts.apply(
            lambda row: calculate_normal_counts(row, controlbam, fasta, window=distance), axis=1
        ).apply(pd.Series)
    else:
        indelcounts['control_alt_counts'] = 0
        indelcounts['control_total_counts'] = 0
        indelcounts['Positions'] = pd.NA
        indelcounts['Distance'] = pd.NA


    # apply filters
    indelcounts = indelcounts[(indelcounts['counts'] >= minreads) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['control_alt_counts'] <= maxcontrol) | (indelcounts['ref']=='.')]
    indelcounts = indelcounts[(indelcounts['Distance'] <= distance) | (indelcounts['ref']=='.')]
    return(indelcounts)

# Conversion function for CSV to hotspot format
def convert_to_hotspot_format(input_csv, output_file):
    with open(input_csv, 'r') as csv_file, open(output_file, 'w', newline='') as out_file:
        reader = csv.DictReader(csv_file)
        writer = csv.writer(out_file)
        
        for row in reader:
            dna_sequence = row['DNA Sequence']
            pam = row['PAM']
            mismatch = row['Mismatch']
            bulge_type = row['Bulge Type']
            bulge_size = row['Bulge Size']
            chromosome = row['Chromosome']
            strand = row['Strand Direction']
            start = row['Start']

            mismatch = mismatch if mismatch else "N/A"
            bulge_type = bulge_type if bulge_type else ""
            bulge_size = bulge_size if bulge_size else ""
            chrom_pos = f"{chromosome}:{strand}{start}"

            writer.writerow([dna_sequence, pam, '', mismatch, '', chrom_pos])

# Convert CSV to hotspot format if it has the specified header
def convert_csv_if_needed(input_file):
    with open(input_file, 'r') as file:
        first_line = file.readline().strip()
        expected_columns = "Source,DNA Sequence,PAM,Chromosome,Strand Direction,Start,Bulge Type,Mismatch,Bulge Size,On_target"
        if first_line == expected_columns:
            output_file = input_file + '.converted'
            convert_to_hotspot_format(input_file, output_file)
            return output_file
    return input_file

def main():

    parser = argparse.ArgumentParser(description='Find indels in a bam file at BED coordinates')
    parser.add_argument('-f','--fasta',type=str,default="/storage1/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh_v4.3.6/hg38_PLVM_CD19_CARv4_cd34.fa",help='Reference fasta file')
    parser.add_argument('-w','--window',type=int,default=100,help='Distance between off-target sites for merging intervals.')
    parser.add_argument('-d','--distance',type=int,default=25,help='Window size around off-target site to identify mutations')
    parser.add_argument('-m','--minreads',type=int,default=1,help='Minimum supporting reads to report an indel/bnd event.')
    parser.add_argument('-c','--maxcontrol',type=int,default=0,help='Maximum supporting reads to to report an indel/bnd event.')
    parser.add_argument('-s','--strict',action='store_true',help='Exit if a read cannot be properly parsed.')
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

    # open bam file(s)
    expsamfile = pysam.AlignmentFile(args.expbamfile,"rc",reference_filename=args.fasta)
    consamfile = pysam.AlignmentFile(args.conbamfile,"rc",reference_filename=args.fasta)
    # open fasta file
    refFasta = pysam.FastaFile(args.fasta)

    # print to outfile or stdout
    if args.outfile:
        sys.stdout = open(args.outfile, 'w')

    print("\t".join('chrom start end pam_positions total_reads indel_reads indel_fraction control_reads control_indel_reads control_indel_fraction indel_count indel_info bnd_count bnd_info target_info is_target'.split(' ')))

    # iterate over mergedBedPr intervals:
    for index, row in mergedBedDf.iterrows():

        indels = get_indels(bam=expsamfile,controlbam=consamfile,chr=row['Chromosome'],start=row['Start'],end=row['End'],
                            fasta=refFasta,window=args.window,distance=args.distance,pam_positions=row['Pos'],minreads=args.minreads,maxcontrol=args.maxcontrol)
        
        total_reads, indel_reads, control_total_reads, control_indel_reads = 0, 0, 0, 0
        indel_fraction, control_indel_fraction = 0, 0
        indel_keys, bnd_keys = '.', '.'
        bnds = pd.DataFrame()
        offtargetsites = ';'.join(row['Info']) if len(row['Info']) > 0 else '.'
        ontarget = row['Ontarget']


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

        print(f"{row['Chromosome']}\t{row['Start']}\t{row['End']}\t{';'.join([ str(x) for x in row['Pos'] ])}\t{total_reads}\t{indel_reads}\t{indel_fraction}\t{control_total_reads}\t{control_indel_reads}\t{control_indel_fraction}\t{len(indels)}\t{indel_keys}\t{len(bnds)}\t{bnd_keys}\t{offtargetsites}\t{ontarget}")



    # close stdout or outfile
    if args.outfile:
        sys.stdout.close()

    expsamfile.close()
    consamfile.close()
    refFasta.close()


if __name__ == "__main__":
    main()
