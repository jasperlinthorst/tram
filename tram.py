#!/usr/bin/env python

import pysam
import os
import gzip
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

def rightinformativeclipping(read,trfstart):
    for rqp,rrp in read.get_aligned_pairs(matches_only=True):
        if rrp>=trfstart:
            return read.infer_query_length()-rqp

def leftinformativeclipping(read,trfend):
    for rqp,rrp in read.get_aligned_pairs(matches_only=True):
        if rrp>=trfend:
            return rqp

def estimateexpansion(pysamfile,chrom,trfstart,trfend,wiggle=500,returnseq=False,usesplitreads=True):
    #extract reads with more than one alignment
    mreads=dict()
    for read in pysamfile.fetch(chrom,trfstart,trfend):
        if (read.reference_start<trfstart-wiggle and read.reference_end>trfstart) or \
            (read.reference_start<trfend and read.reference_end>trfend+wiggle):
            if read.qname in mreads:
                mreads[read.qname].append(read)
            else:
                mreads[read.qname]=[read]
    se=[]
    pe=[]
    lce=[]
    rce=[]

    salignments=[]
    palignments=[]
    calignments=[]

    sorientations=[]
    porientations=[]
    lcorientations=[]
    rcorientations=[]
    
    sseq=[]
    pseq=[]
    rcseq=[]
    lcseq=[]
    
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
            if reads[0].is_supplementary or reads[0].is_secondary:
                continue
            
            if reads[0].reference_start<trfstart-wiggle and reads[0].reference_end>trfend+wiggle:
                qa,qb=querypositionsfromsingleread(reads[0],trfstart,trfend)
                se.append(qb-qa)
                sorientations.append(reads[0].is_reverse)
                if returnseq:
                    sseq.append(s[qa:qb])
                salignments.append(reads[0])
            
            else:
                if reads[0].reference_start<trfstart-wiggle and ((reads[0].cigar[-1][0]==4 or reads[0].cigar[-1][0]==5) and reads[0].cigar[-1][1]+reads[0].reference_end>trfend):
                    ric=rightinformativeclipping(reads[0],trfstart)
                    rce.append(ric)
                    rcorientations.append(reads[0].is_reverse)
                    if returnseq:
                        rcseq.append(s[-ric:])
                    calignments.append(reads[0])

                if reads[0].reference_end>trfend+wiggle and ((reads[0].cigar[0][0]==4 or reads[0].cigar[0][0]==5) and reads[0].reference_start-reads[0].cigar[0][1]<trfstart):
                    lic=leftinformativeclipping(reads[0],trfend)
                    lce.append(lic)
                    lcorientations.append(reads[0].is_reverse)
                    if returnseq:
                        lcseq.append(s[:lic])
                    calignments.append(reads[0])
                    
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
            pe.append(qb-qa)
            porientations.append(read1.is_reverse)
            palignments.append((read1,read2))
            if returnseq:
                pseq.append(s[qa:qb])
    
    #if len(se)>100: #exclude estimates that are based on more then this many reads
    #    print "Coverage too high for proper estimate",len(e)
    #    return None
    
    return se,sorientations,salignments,sseq,pe,porientations,palignments,pseq,lce,lcorientations,lcseq,rce,rcorientations,rcseq,calignments

def main():
    desc="""
    Type 'tram <positional_argument> --help' for help on a specific subcommand.
    """
    
    parser = argparse.ArgumentParser(prog="tram", usage="tram -h for usage", description=desc)
    parser.add_argument('bamfile', type=str, help='BAM file that contains the long read alignments.')
    parser.add_argument('bedfile', type=str, help='BED file specifying the (repeat) regions to scan for expansion/contraction.')
    parser.add_argument("--wiggle", dest="wiggle", type=int, default=None, help="How far split reads may span region start and end point and still be considered for length estimation (default: equal to length of repeat pattern).")
    parser.add_argument("--nosplitreads", dest="usesplitreads", action="store_false", default=True, help="Don't use splitreads for length estimation.")
    parser.add_argument("--seq", dest="returnseq", action="store_true", default=False, help="Return the sequence in between the two flanks.")
    parser.add_argument("--slice", dest="bamslice", action="store_true", default=False, help="Produce bam files (s=singleread,p=splitreads,c=clipping) that contain subsets of the alignments that were used for estimating the lengths.")
    
    args = parser.parse_args()
    
    bamfile=pysam.AlignmentFile(args.bamfile, "rb")
    
    if not args.bamfile.endswith(".bam"):
        print "Invalid bam file."
        return
    
    if args.bamslice:
        sslicefile=pysam.AlignmentFile(os.path.basename(args.bamfile).replace(".bam",".sslice.bam"), "wb", template=bamfile)
        pslicefile=pysam.AlignmentFile(os.path.basename(args.bamfile).replace(".bam",".pslice.bam"), "wb", template=bamfile)
        cslicefile=pysam.AlignmentFile(os.path.basename(args.bamfile).replace(".bam",".cslice.bam"), "wb", template=bamfile)
    
    if args.wiggle!=None:
         w=args.wiggle
    
    try:
        with open(args.bedfile) as trf:
            for i,line in enumerate(trf):
                cols=line.split("\t")
                chrom,trfstart,trfend=cols[0],int(cols[1]),int(cols[2])
                
                if args.wiggle==None:
                    w=trfend-trfstart

                v=estimateexpansion(bamfile,chrom,trfstart,trfend,wiggle=w,usesplitreads=args.usesplitreads,returnseq=args.returnseq)
                
                se,sorientations,salignments,sseq,pe,porientations,palignments,pseq,lce,lcorientations,lcseq,rce,rcorientations,rcseq,calignments=v
                
                sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s"%(chrom,trfstart,trfend,trfend-trfstart,\
                                                            ",".join([str(x) for x in se]),\
                                                            ",".join([str(x) for x in pe]),\
                                                            ",".join([str(x) for x in lce]),\
                                                            ",".join([str(x) for x in rce]),\
                                                            ))
                
                if args.returnseq:
                    sys.stdout.write("\t%s"%",".join(sseq))
                    sys.stdout.write("\t%s"%",".join(pseq))
                    sys.stdout.write("\t%s"%",".join(lcseq))
                    sys.stdout.write("\t%s"%",".join(rcseq))
                
                sys.stdout.write("\n")
                
                if args.bamslice:
                    for a in salignments:
                        if isinstance(a,tuple):
                            sslicefile.write(a[0])
                            sslicefile.write(a[1])
                        else:
                            sslicefile.write(a)
                    for a in palignments:
                        if isinstance(a,tuple):
                            pslicefile.write(a[0])
                            pslicefile.write(a[1])
                        else:
                            pslicefile.write(a)
                    for a in calignments: 
                        if isinstance(a,tuple):
                            cslicefile.write(a[0])
                            cslicefile.write(a[1])
                        else:
                            cslicefile.write(a)
    
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass
    
    if args.bamslice:
        sslicefile.close()
        pysam.sort("-o", os.path.basename(args.bamfile).replace(".bam",".sslice.bam"), os.path.basename(args.bamfile).replace(".bam",".sslice.bam"))
        pysam.index(os.path.basename(args.bamfile).replace(".bam",".sslice.bam"))

        pslicefile.close()
        pysam.sort("-o", os.path.basename(args.bamfile).replace(".bam",".pslice.bam"), os.path.basename(args.bamfile).replace(".bam",".pslice.bam"))
        pysam.index(os.path.basename(args.bamfile).replace(".bam",".pslice.bam"))
        
        cslicefile.close()
        pysam.sort("-o", os.path.basename(args.bamfile).replace(".bam",".cslice.bam"), os.path.basename(args.bamfile).replace(".bam",".cslice.bam"))
        pysam.index(os.path.basename(args.bamfile).replace(".bam",".cslice.bam"))
        

    bamfile.close()
