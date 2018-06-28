#!/usr/bin/env python

import pysam
import os
import gzip
import random
from intervaltree import IntervalTree, Interval
import argparse
import sys, errno
import numpy as np

from scipy.optimize import minimize
from scipy.cluster.vq import kmeans2
from scipy.stats import norm

import warnings

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

def bimodal_ll(params,data):
    mu1,mu2,s=params
    return -np.sum(np.log(np.amax(np.column_stack((norm.pdf(data, loc=mu1, scale=s), norm.pdf(data, loc=mu2, scale=s))),axis=1)+.00000001))

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
    parser.add_argument("--plot", dest="plot", action="store_true", default=False, help="Create a plot for every locus in the bedfile.")
    parser.add_argument("-i", dest="interactive", action="store_true", default=False, help="Show interactive plot, (will pause until window is closed).")
    
    args = parser.parse_args()
    
    bamfile=pysam.AlignmentFile(args.bamfile, "rb")
    
    prefix=os.path.basename(args.bamfile).replace(".bam","")

    if not args.bamfile.endswith(".bam"):
        print "Invalid bam file."
        return
    
    if args.bamslice:
        slicefile=pysam.AlignmentFile(prefix+".slice.bam", "wb", template=bamfile)
        sslicefile=pysam.AlignmentFile(prefix+".sslice.bam", "wb", template=bamfile)
        pslicefile=pysam.AlignmentFile(prefix+".pslice.bam", "wb", template=bamfile)
        cslicefile=pysam.AlignmentFile(prefix+".cslice.bam", "wb", template=bamfile)
    
    if args.wiggle!=None:
         w=args.wiggle
    
    if args.plot:
        from matplotlib import pyplot as plt

    try:
        slicedalignments=set()
        sslicedalignments=set()
        pslicedalignments=set()
        cslicedalignments=set()

        with open(args.bedfile) as trf:
            #write a header that describes the columns
            cols=['chrom','trfstart','trfend','reflength','length estimate 1','length estimate 2','within read alignment estimates','split read alignment estimates','left clipping','right clipping']
            if args.returnseq:
                cols+=['within read alignment estimates (sequence)','split read alignment estimates (sequence)','left clipping (sequence)','right clipping (sequence)']

            sys.stdout.write("#"+"\t".join(cols)+"\n")
            for i,line in enumerate(trf):
                cols=line.rstrip().split("\t")
                chrom,trfstart,trfend=cols[0],int(cols[1]),int(cols[2])
                
                if len(cols)>=4:
                    name=cols[3]
                else:
                    name="unknown"
                
                if args.wiggle==None:
                    w=trfend-trfstart
                
                v=estimateexpansion(bamfile,chrom,trfstart,trfend,wiggle=w,usesplitreads=args.usesplitreads,returnseq=args.returnseq)
                
                se,sorientations,salignments,sseq,pe,porientations,palignments,pseq,lce,lcorientations,lcseq,rce,rcorientations,rcseq,calignments=v
                
                #Fit bimodal gaussion distribution to come up with diploid assignment
                lengthdist=np.array(se+pe).astype(float)
                
                fit=None
                if len(lengthdist)>1:
                    with warnings.catch_warnings(): #to prevent scipy kmeans userwarnings
                        warnings.simplefilter("ignore")
                        s=set([1])
                        while len(s)==1:
                            initmu,l=kmeans2(lengthdist,2,iter=10) #use kmeans to initialize parameters
                            s=set(l)
                    fit=minimize(bimodal_ll,x0=(initmu[0],initmu[1],1),bounds=[(None,None),(None,None),(1,None)],method='L-BFGS-B',args=lengthdist) #minimize negative likelihood for bimodal gaussian distribution
                    lengthestimates=sorted([int(fit.x[0]),int(fit.x[1])])
                elif len(lengthdist)==1:
                    lengthestimates=(lengthdist[0],lengthdist[0])
                else:
                    lengthestimates=(None,None)
                
                # ll=bimodal_ll(fit.x, lengthdist)*-1
                
                sys.stdout.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s"%(chrom,trfstart,trfend,trfend-trfstart,lengthestimates[0],lengthestimates[1], \
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
                
                if args.plot:
                    plt.clf()
                    # plt.subplot(211)
                    # plt.plot([1]*len(lengthdist),lengthdist,'*')
                    # plt.axhline(y=fit.x[0],linewidth=1,color='black',linestyle='solid')
                    # plt.axhline(y=fit.x[1],linewidth=1,color='black',linestyle='solid')
                    plt.title(" ".join(cols[:3]))
                    # plt.subplot(212)
                    #plt.hist(lengthdist, bins=25, density=True, alpha=0.6, color='g')
                    plt.hist(lengthdist, bins=25, density=True, alpha=0.6, color='g')
                    xmin, xmax = plt.xlim()
                    x = np.linspace(xmin, xmax, 100)
                    
                    if fit!=None:
                        p1 = norm.pdf(x, lengthestimates[0], fit.x[2])
                        plt.plot(x, p1, 'k', linewidth=2)
                        p2 = norm.pdf(x, lengthestimates[1], fit.x[2])
                        plt.plot(x, p2, 'k', linewidth=2)
                    
                    if lengthestimates[0]!=None:
                        plt.axvline(x=lengthestimates[0],linewidth=1,color='b',linestyle='-')
                    
                    if lengthestimates[1]!=None:
                        plt.axvline(x=lengthestimates[1],linewidth=1,color='r',linestyle='-')
                    
                    if args.interactive:
                        plt.show()
                    else:
                        plt.savefig(prefix+"_"+"_".join(cols[:3])+"_"+name+".png")

                if args.bamslice:
                    
                    for a in salignments:
                        aid=(a.query_name,a.pos)
                        if aid not in sslicedalignments:
                            sslicefile.write(a)
                            sslicedalignments.add(aid)
                            
                            if aid not in slicedalignments:
                                slicefile.write(a)
                                slicedalignments.add(aid)

                    for pa in palignments:
                        for a in pa:
                            aid=(a.query_name,a.pos)
                            if aid not in pslicedalignments:
                                pslicefile.write(a)
                                pslicedalignments.add(aid)
                                
                                if aid not in slicedalignments:
                                    slicefile.write(a)
                                    slicedalignments.add(aid)

                    for a in calignments:
                        aid=(a.query_name,a.pos)
                        if aid not in cslicedalignments:
                            cslicefile.write(a)
                            cslicedalignments.add(aid)
                        
                        if aid not in slicedalignments:
                            slicefile.write(a)
                            slicedalignments.add(aid)
    
    except IOError as e:
        if e.errno == errno.EPIPE:
            pass
    
    if args.bamslice:
        
        slicefile.close()
        pysam.sort("-o", os.path.basename(args.bamfile).replace(".bam",".slice.bam"), os.path.basename(args.bamfile).replace(".bam",".slice.bam"))
        pysam.index(os.path.basename(args.bamfile).replace(".bam",".slice.bam"))
        
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
