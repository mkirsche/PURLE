# PURLE

Here we introduce PURLE, software for Probabilistic Uncontained Read Length Estimation.  
  
Given a set of genomic reads, especially in higher coverage settings, it is expected that many of them will be strictly contained in other reads (i.e., the interval of the genome they come from is a subset of the interval of some other read), and hence give mostly redundant information.  Given a length distribution of reads from a genome with a given length, this software 1) Gets the expected distribution, or optionally a sampled distribution, of read lengths for the uncontained reads (assuming each read independently comes from a uniformly random position in the genome) and outputs a file containing these read lengths and 2) Helps in an estimating a length cutoff by calculating the expected number of uncontained reads below a predetermined set of length thresholds.

## Compilation: 

javac LengthDistribution.java

## Usage: 

java LengthDistribution Inputfile Outputfile Genomelength [--sample]  
&nbsp;&nbsp;&nbsp;&nbsp;Inputfile should contain the length of one sequence (in bp) per line  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Alternately, the input file can be in FASTQ or FASTA format.  
&nbsp;&nbsp;&nbsp;&nbsp;Outputfile is the file to output sample noncontained read lengths to  
&nbsp;&nbsp;&nbsp;&nbsp;Genomelength is an integer representing the length (in bp) of the genome  
&nbsp;&nbsp;&nbsp;&nbsp;--sample is an optional flag to produce a sample length distribution instead of the average




