# TRAM (Tandem Repeat Assessment Method)

TRAM is a simple utility script that measures a basepair distance between two anchor points on a reference sequence for each read spanning both points.
It can be used to assess wheter a repeat element (defined on the reference genome) might be expanded or contracted.
It considers split read pairs in case an aligner breaks on the expansion/contraction.

## INSTALL

**python setup.py install**

Tram depends on the pysam and intervaltree packages.

## USAGE

tram &lt;bam-file&gt; &lt;bed-file&gt;

The bam-file should contain the reads aligned to your reference sequence.
The bed-file should contain three columns specifying the chromosome, the start position and the end position of the repeat element which you'd like to investigate.

## OUTPUT
Output is a tab-separated value list that contains the following columns:

The first columns are essentialy a repetition of what was specified in the bed file.

1. Reference chromosome.
2. Start of the repeat elements
3. End of the repeat elements
4. Length on the reference (essentially 3-2)

The second part contain the (comma separated) length measurements from the sequence reads:

5. Length measurements of alignments that span the locus
6. Length measurements of split-read alignments that break on the flanks of the specified region
7. Length of the sequence that was clipped on the left flank, for alignments that did not span the locus
8. Length of the sequence that was clipped on the rihgt flank, for alignments that did not span the locus

When --seq was specified the following additional columns are output:

9. Sequence that corresponds to the length measurements in column 5
10. Sequence that corresponds to the length measurements in column 6
11. Sequence that corresponds to the length measurements in column 7
12. Sequence that corresponds to the length measurements in column 8

For more info, run "tram -h":

>Type 'tram <positional_argument> --help' for help on a specific subcommand.
>positional arguments:
> * bamfile          BAM file that contains the long read alignments.
> * bedfile          BED file specifying the (repeat) regions to scan for expansion/contraction.
>optional arguments:
> * -h, --help       show this help message and exit
> * --wiggle WIGGLE  How far split reads may span region start and end point and still be considered for length estimation (default: equal to length of repeat pattern).
> * --nosplitreads   Don't use splitreads for length estimation.
> * --seq            Return the sequence in between the two flanks.
> * --slice          Produce bam files (s=singleread,p=splitreads,c=clipping) that contain subsets of the alignments that were used for estimating the lengths.