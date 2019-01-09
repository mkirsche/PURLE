# PURLE

Here we introduce PURLE, software for Probabilistic Uncontained Read Length Estimation.

## Compilation: 

javac LengthDistribution.java

## Usage: 

java LengthDistribution Inputfile Outputfile Genomelength
  Inputfile should contain the length of one sequence (in bp) per line
   Alternately, the input file can be in FASTQ or FASTA format.
  Outputfile is the file to output sample noncontained read lengths to
  Genomelength is an integer representing the length (in bp) of the genome




