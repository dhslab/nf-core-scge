#!/usr/bin/python3

from __future__ import division
import pysam
import argparse
import re
import string
import csv
import scipy.stats as stats

# reverse complement function
def revcomp(seq):
    tab = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh') # maketrans <- maps the reverse complement
    return seq.translate(tab)[::-1] # translate(x)[::-1] <- works backward through the string, effectively reversing the string

# function to get indels for one amplicon from bam file
def get_indels(bam,chr,start,end,pampos,strand,srchSeq,srchSeqMM):

    dat = {'all':0,'indels':{},'sv':{},'seq':0}
    
    for read in bam.fetch(chr, start, end, multiple_iterators = True):

        # skip if not primary alignment or a duplicate or alignment doesnt overlap start,end
        if read.is_mapped is False or \
            read.is_duplicate is True or \
            read.is_secondary is True or \
            read.is_supplementary is True or \
            read.reference_start > end or \
            read.reference_end < start:
            continue

        if read.mapping_quality == 0:
            continue

        if read.has_tag('NM') is True and int(read.get_tag('NM'))>5:
            continue

        # check for a chimeric read
        if read.is_proper_pair is False and read.mate_is_mapped is True:
            mate_contig = read.next_reference_name
            mate_start = read.next_reference_start
            if mate_contig != read.reference_name or abs(mate_start - read.reference_start) > 1000:
                dat['sv'][':'.join([read.reference_name,str(read.reference_start),read.next_reference_name,str(read.next_reference_start)])] = 1
            continue

        if read.has_tag('SA') is True:
            sChr, sPos, sStrand, sCigar, sMq, sNm = read.get_tag('SA').split(',')
            if sChr != read.reference_name and int(sMq)>0:
                dat['sv'][':'.join([read.reference_name,str(read.reference_start),sChr,str(sPos)])] = 1
            continue

        # check that read start and end are within w bp of the amplicon target start and end (to avoid spurious alignments--skip if this is not true) 
        cigar = read.cigartuples # get cigar info

        leftSoftClip = 0
        rightSoftClip = 0
        if cigar[0][0] == 4:
            leftSoftClip = cigar[0][1]

        if cigar[-1][0] == 4:
            rightSoftClip = cigar[-1][1]

        # skip if read is not aligned across the pampos (incl. soft clips)
#        if read.is_proper_pair is True and (read.reference_start - leftSoftClip > pampos or read.reference_end + rightSoftClip < pampos):
#            continue

        # skip if read doesnt fully span the pampos
        if read.reference_start >= pampos or read.reference_end <= pampos:
            continue

        # we cant resolve these so will skip for now.
        if leftSoftClip > 0 or rightSoftClip > 0:
            continue

        # count all reads covering this bed record position
        dat['all'] += 1

        # remove start and end cigar operations that are N/H/S/X (ie, not M, D, or I).
        # These are soft or hard clips and not aligned portions of the read
        if cigar[0][0] > 2:
            del cigar[0]
            
        if cigar[-1][0] > 2:
            del cigar[-1]                    
                
        if len(cigar) == 1 and cigar[0][0] == 0: # no indels, so continue
            continue                       

        # get reference position start for the alignment
        # then skip beginning and trailing aligned blocks that are aligned
        # to the reference
        rpos = read.reference_start

        var = [] # operation: [I]nsertion, [D]eletion, or [C]omplex
        varlen = [] # length of indel
        vartype = [] # variant type
        varst = [] # ref pos of variant start 
        varen = [] # ref pos of variant end
        deltabases = 0 # sum of bases inserted/deleted

        for cg in cigar: # rebuild cigar and keep track
            s = rpos
            e = rpos
            if cg[0] == 0 or cg[0] == 2 or cg[0] == 3:
                e = rpos + cg[1]                
            elif cg[0] == 1:
                e = rpos + 1
                
            # do not include the cigar operation if it does not overlap start,end...
            if e < start or s > end:
                rpos = e
                continue

            # ...or its a matching segment that goes beyond start,end since these dont matter            
            elif cg[0] == 0 and (s < start or e > end):
                rpos = e
                continue
            
            elif cg[0] == 0: # aligned block
                var.append(str(cg[1]) + "M")
                varlen.append(cg[1])
                vartype.append("M")
                varst.append(s)
                varen.append(e)
                
            elif cg[0] == 1: # an insertion
                var.append(str(cg[1]) + "I")
                vartype.append("I")
                varlen.append(cg[1])
                varst.append(s)
                varen.append(e)
                deltabases += cg[1]
               
            elif cg[0] == 2 or cg[0] == 3: # a deletion or reference skip
                var.append(str(cg[1]) + "D")
                vartype.append("D")
                varlen.append(cg[1])
                varst.append(s)
                varen.append(e)
                deltabases += cg[1]
                
            else: # account for other operations
                var.append(str(cg[1]) + "X")
                vartype.append("C")
                varlen.append(cg[1])
                varst.append(s)
                varen.append(e)
                deltabases += cg[1]
                
            rpos = e # increment reference position

        # do not continue if there are no indels
        if len(var) == 0:
            continue
            
        # combine multiple indels into a complex one
        if vartype[0] == 'M':
            var = var[1:]
            varlen = varlen[1:]
            vartype = vartype[1:]
            varst = varst[1:]
            varen = varen[1:]

        if vartype[-1] == 'M':
            var = var[:-1]
            varlen = varlen[:-1]
            vartype = vartype[:-1]
            varst = varst[:-1]
            varen = varen[:-1]
        
        if len(var) > 1:
            var = ["".join(var)]
            varlen = [sum(varlen)]
            vartype = ["C"]
            varst = [min(varst)]
            varen = [max(varen)]
            
        i = 0
        # get position relative to PAM, strand
        relpos = None
        if (vartype[i] == "D" or vartype[i] == "C") and varst[i] <= pampos and varen[i] >= pampos:
            relpos = 0
                
        elif strand == "+" and varst[i] >= pampos: # fwd downstream
            relpos = varst[i] - pampos
                
        elif strand == "+" and varen[i] <= pampos: # fwd upstream
            relpos = varen[i] - pampos
    
        elif strand == "-" and varen[i] <= pampos: # rev downstream
            relpos = pampos - varen[i]
    
        elif strand == "-" and varst[i] >= pampos: # rev upstream
            relpos = pampos - varst[i] 

        # key is reference_position:cigar:vartype:len:relpos
        key = str(varst[i]) + ":" + var[i] + ":" + vartype[i] + ":" + str(varlen[i]) + ":" + str(relpos) 
            
        if key not in dat['indels'].keys():
            dat['indels'][key] = 0

        dat['indels'][key] += 1 # count the reads with this variant        

    return(dat)


srchSeq = 'TTAATTGAGTTGTCATATGTTAATAACGGT' # this is the dsODN sequence
srchSeqMM = 2 # mismatches tolerated when looking for the search sequence
maxInControl = 10

parser = argparse.ArgumentParser(description='Find indels in a bam file at BED coordinates')
parser.add_argument('bed',help='Amplicon bed file')
parser.add_argument('expbamfile',help='BAM file')
parser.add_argument('conbamfile',help='BAM file')

args = parser.parse_args()

bed = args.bed

# open bam file(s)
expsamfile = pysam.AlignmentFile(args.expbamfile,"rc",reference_filename="/storage1/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh/hg38_PLVM_CD19_CARv4_cd34.fa")
consamfile = pysam.AlignmentFile(args.conbamfile,"rc",reference_filename="/storage1/fs1/dspencer/Active/clinseq/projects/scge/data/refdata/singh/hg38_PLVM_CD19_CARv4_cd34.fa")

distup = 20
distdown = 20

print("\t".join('id chrom pos strand edited_reads edited_indels edited_unique_indels control_reads control_indels control_unique_indels indel_info exp_svs con_svs'.split(' ')))
#rec = [contig,pampos,strand,datExp['all'],expindels,sum(datExp['filtindels'].values()),datExp['seq'],datCon['all'],controlindels,sum(datCon['filtindels'].values()),datCon['seq'],";".join(map(str,vsfilt)) or '.']
with open(bed) as f:
    for line in f:
        rec = line.strip().split(',')

        p = re.compile("^(\S+)\:([+-])(\d+)")
        m = p.search(rec[5])
        contig = m[1]
        pampos = int(m[3])
        strand = m[2]
        if strand == '+':
            pampos = pampos + 20

        st = pampos - distup
        en = pampos + distdown
        if strand == "-":
            st = pampos - distdown
            en = pampos + distup

        datExp = get_indels(expsamfile,contig,st,en,pampos,strand,srchSeq,srchSeqMM)
        datCon = get_indels(consamfile,contig,st,en,pampos,strand,srchSeq,srchSeqMM)

        datExp['filtindels'] = {}
        datExp['filtseq'] = {}
        datCon['filtindels'] = {}
        datCon['filtseq'] = {}

        expindels = 0
        controlindels = 0
        
        for k in set(list(datExp['indels'].keys()) + list(datCon['indels'].keys())): # iterate through all exp and con indels

            if (datCon['indels'].get(k) or 0) > maxInControl: # do not report indels with >maxInControl reads--these are inherited or artifacts. delete from dicts
                if datExp['indels'].get(k):
                    del datExp['indels'][k]
                    
                if datCon['indels'].get(k):
                    del datCon['indels'][k]

            else:
                datExp['filtindels'][k] = datExp['indels'].get(k) or 0
                if datExp['filtindels'][k] > 0:
                    expindels += 1
                    
                datCon['filtindels'][k] = datCon['indels'].get(k) or 0
                if datCon['filtindels'][k] > 0:
                    controlindels += 1


        # assemble filtered data
        vsfilt = []
        for k in datExp['filtindels'].keys():
            vsfilt.append(":".join(map(str,[k,datExp['filtindels'].get(k) or 0,datCon['indels'].get(k) or 0])))

        expSvs = ';'.join(datExp['sv'].keys()) or '.'
        controlSvs = ';'.join(datCon['sv'].keys()) or '.'
        rec = [rec[5],contig,pampos,strand,datExp['all'],expindels,sum(datExp['filtindels'].values()),datCon['all'],controlindels,sum(datCon['filtindels'].values()),";".join(map(str,vsfilt)) or '.', expSvs, controlSvs]
        
        print("\t".join(map(str,rec)))

f.close()
expsamfile.close()
consamfile.close()
