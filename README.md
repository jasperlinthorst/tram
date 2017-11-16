# TRAM (Tandem Repeat Assessment Method)

TRAM is a simple utility script that measures a basepair distance between two anchor points on a reference sequence for each read spanning both points.
It can be used to assess wheter a repeat element (defined on the reference genome) might be expanded or contracted.
It considers split read pairs in case an aligner break on the expansion/contraction.

## INSTALL

python setup.py install

Tram depends on the pysam and intervaltree packages.

## USAGE

tram <bam-file> <bed-file>

The bam-file should contain the reads aligned to your reference sequence.
The bed-file should contain three columns specifying the chromosome, the start position and the end position of the repeat element which you'd like to ivnestigate.

Output is a tab-separated value list that contains the various measurements.

For more info, run "tram -h".
