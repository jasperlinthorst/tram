#!/usr/bin/env python

import pysam
from matplotlib import pyplot as plt
import os
import gzip
import numpy as np
import random
from intervaltree import IntervalTree, Interval
import argparse
import sys, errno


def querypositionsfromsplitreads(read1,read2,trfstart,trfend):
    r1lhclip=read1.cigar[0][1] if read1.cigar[0][0]==5 else 0
    r2lhclip=read2.cigar[0][1] if read2.cigar[0][0]==5 else 0
    for r1qp,r1rp in read1.get_aligned_pairs(matches_only=True):
        if r1rp>=trfstart:
            break
        else:
            qa=r1qp+r1lhclip
            ra=r1rp
    for r2qp,r2rp in read2.get_aligned_pairs(matches_only=True):
        if r2rp>=trfend:
            break
        else:
            qb=r2qp+r2lhclip
            rb=r2rp
    return qa,qb

def querypositionsfromsingleread(read,trfstart,trfend):
    rlhclip=read.cigar[0][1] if read.cigar[0][0]==5 else 0
    for rqp,rrp in read.get_aligned_pairs(matches_only=True):
        if rrp>=trfstart:
            break
        else:
            qa=rqp+rlhclip
            ra=rrp
    for rqp,rrp in read.get_aligned_pairs(matches_only=True):
        if rrp>=trfend:
            break
        else:
            qb=rqp+rlhclip
            rb=rrp
    return qa,qb

def estimateexpansion(pysamfile,chrom,trfstart,trfend,wiggle=1000,returnseq=False,usesplitreads=True):
    #extract reads with more than one alignment
    mreads=dict()
    for read in pysamfile.fetch(chrom,trfstart,trfend):
        if (read.reference_start<trfstart-wiggle and read.reference_end>trfstart) or \
            (read.reference_start<trfend and read.reference_end>trfend+wiggle):
            if read.qname in mreads:
                mreads[read.qname].append(read)
            else:
                mreads[read.qname]=[read]
    e=[]
    alignments=[]
    usedreads=[]
    orientations=[]
    seq=[]
    #select left and right most alignment that does overlap the trf
    for readname in mreads:
        reads=sorted(mreads[readname],key=lambda r: r.reference_start+r.alen)
        
        if returnseq:
            s=""
            for read in reads:
                if not read.is_supplementary and not read.is_secondary:
                    s=read.query_sequence
                    break
            assert(s!=None)
        
        if len(reads)==1: #estimate length within a single aligned segment
            if not (reads[0].reference_start<trfstart-wiggle and reads[0].reference_end>trfend+wiggle):
                continue
            if reads[0].is_supplementary or reads[0].is_secondary:
                continue
            qa,qb=querypositionsfromsingleread(reads[0],trfstart,trfend)
            e.append(qb-qa)
            alignments.append(reads[0])
            usedreads.append(reads[0].qname)
            orientations.append(reads[0].is_reverse)
        else: #handle split-read alignments, take left and rightmost aligned segment wrt the query
            if not usesplitreads:
                continue
            read1=reads[0]
            read2=reads[-1]
            if read1.reference_start>trfstart-wiggle or read1.reference_start+read1.alen>trfend+wiggle or \
                read2.reference_start<trfstart-wiggle or read2.reference_start+read2.alen<trfend+wiggle:
                continue
            if not read1.is_reverse==read2.is_reverse:
                continue
            qa,qb=querypositionsfromsplitreads(read1,read2,trfstart,trfend)
            e.append(qb-qa)
            alignments.append((read1,read2))
            usedreads.append(read1.qname)
            orientations.append(read1.is_reverse)
        
        if returnseq:
            seq.append(s[qa:qb])
    
    if len(e)==0:
        return None
    if len(e)>100: #exclude estimates that are based on more then this many reads
        print "Coverage too high for proper estimate",len(e)
        return None
    
    qlength=int(np.median(e))
    
    return qlength,usedreads,alignments,e,seq

def main():
    desc="""
    Type 'tram <positional_argument> --help' for help on a specific subcommand.
    """
    
    parser = argparse.ArgumentParser(prog="tram", usage="tram -h for usage", description=desc)
    parser.add_argument('bamfile', type=str, help='BAM file that contains the long read alignments.')
    parser.add_argument('bedfile', type=str, help='BED file specifying the (repeat) regions to scan for expansion/contraction.')
    parser.add_argument("--wiggle", dest="wiggle", type=int, default=500, help="How far split reads may span region start and end point and still be considered for length estimation.")
    parser.add_argument("--nosplitreads", dest="usesplitreads", action="store_false", default=True, help="Don't use splitreads for length estimation.")
    parser.add_argument("--slice", dest="slice", action="store_true", default=False, help="Produce a bam file that contains only the subset of the alignments that were used for the length estimations.")
    
    args = parser.parse_args()

    bamfile=pysam.AlignmentFile(args.bamfile, "rb")
    
    if not args.bamfile.endswith(".bam"):
        print "Invalid bam file."
        return
    
    if slice:
        slicefile=pysam.AlignmentFile(os.path.basename(args.bamfile).replace(".bam",".slice.bam"), "wb", template=bamfile)
    
    try:
        with open(args.bedfile) as trf:
            for i,line in enumerate(trf):
                cols=line.split("\t")
                chrom,trfstart,trfend=cols[0],int(cols[1]),int(cols[2])
                v=estimateexpansion(bamfile,chrom,trfstart,trfend,wiggle=args.wiggle,usesplitreads=args.usesplitreads)
                
                if v==None:
                    ml,usedreads,alignments,e,seq=(-1,[],[],[],[])
                else:
                    ml,usedreads,alignments,e,seq=v
                
                sys.stdout.write("%s\t%d\t%d\t%d\t%d\t%s\n"%(chrom,trfstart,trfend,trfend-trfstart,ml,",".join([str(x) for x in e])))
                
                if slice:
                    for a in alignments:
                        if isinstance(a,tuple):
                            slicefile.write(a[0])
                            slicefile.write(a[1])
                        else:
                            slicefile.write(a)
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass
    
    if slice:
        slicefile.close()

    bamfile.close()
