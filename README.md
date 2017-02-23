# MConsensus.pl

  Generate a consensus for a specified gene in a specified taxa

## A word of warning

  This script is a PITA. It requires a lot of babysitting of the output to make sure it's outputting something 
  meaningful and not a load of junk. If you're not comfortable looking through the perl code to figure things
  out I would recommend steering clear. If you do try to use it I would recommend at a minimum using the 
  --verbose option and watching the output. Also, look over your alignment output to be sure that it looks 
  like a decent alignment. If your alignment is bad, your consensus is going to be bad. There are a lot of 
  poorly annotated sequences on NCBI and they can mess up everything pretty badly in this script. So, I would 
  recommend not using this script unless you feel like wading into a bit of a mess potentially.
    
## perl MConsensus.pl [options]

  This script takes in an organism name (any scientific taxon) and a gene (or other genetic element) name.
  The script then downloads all the sequences for that query from the NCBI nt database. All sequences are 
  then aligned using mafft. The consensus sequence of the alignment is then output as both a fasta and as 
  a table of counts for each A/T/G/C at each position. 

## Dependencies

### This program is designed to work on Unix 
    
### mafft: mafft is used to align the sequences obtained and needs to be in your $PATH


## Output files

allSeqsAligned.fasta
- All sequences (both blast results and sequences retrieved from NCBI), aligned.

allSeqs.fasta
- All sequences (both blast results and sequences retrieved from NCBI).

Consensus.fasta
- Final file that you will want to look at. Fasta file of the consensus.

ConsensusTable.txt
- Table of the numbers used to generate the consensus. 
- Columns are: 
  * Base - base number 0:length(sequence)
  * A - number of "A"s in the aligned sequence
  * T - number of "T"s in the aligned sequence       
  * G - number of "G"s in the aligned sequence       
  * C - number of "C"s in the aligned sequence
  * Consensus - consensus call - note that the consensus takes the minAF option into account   
  * Missing - number of sequences with zero coverage at this location

organismGis.txt
- List of all GIs matching the query organism. Used to limit the blast search to just the target organism.

originalGis.fasta
- Sequences of the GIs retrieved from NCBI

originalGis.html
- Data retrieved from NCBI search. Used to extract the GIs for sequence retrieval.

originalGis.txt
- GIs extracted from originalGis.html

## Options

--retmax (1000)

  The maximum number of hits to retrieve from the NCBI database matching the organism and gene

--organism ("salamanders")

  The organism to be queried. Use scientific names (e.g. canidae) or other names recognized by 
  NCBI database (e.g. salamanders). I would *highly recommend* testing your search: 
      <http://www.ncbi.nlm.nih.gov/nuccore/>

--gene ("mitochondria+\"complete+genome\"")

  The gene sequence (or other genetic element) to be queried. If you use multiple words, put them 
  in quotes and put a plus between the words. I would *highly recommend* testing your search: 
      <http://www.ncbi.nlm.nih.gov/nuccore/>

--email (matthewvc1@gmail.com)

  Your email. This is used in the query to NCBI.

--outDir (organism_date_time)

  The directory to write the files out to. By default, the directory is named for the organism and
  month/day/time the script is run.

--processors (1)

  The number of processors to use.

--minAF (20)

  This is the minimum allele frequency (in percentage) to include when generating the consensus sequence.
  A lower value will result in more ambiguous bases. 

--allowPartial

  This is a flag that tells the program to allow partial sequences to be used in the alignment and consensus.
  A partial sequence is defined as having either ">" or "<" in either of the gene position fields in the 
  annotation. This option will force the consensus generation algorithm to ignore positions with lots of gaps.
  This means that positions covered by at least one base will be included. Use the consensusTable.txt file 
  to filter out positions with low sequence coverage.

--blackList

  This is a simple list of accessions to exclude from the alignment and consensus. This is found after the
  ">" in the fasta file and up to (but not including) the first space. eg: KY119968.1.

--minLen (0)

  Minimum sequence length to be kept.

--maxLen (INF)

  Maximum sequence length to be kept.

--verbose

  This option prints out status updates while the program runs.

--help 

  Print out this "helpful" information. 


## Examples

You can run it with default settings like this:

`perl MConsensus.pl`

Or you can specify the organism and gene like this:

`perl MConsensus.pl --out cytbCephalopoda --processors 20 --retmax 2000 --organism Cephalopoda --geneName "gene\tcytb,gene\tcob" --verbose`

Or you can tell it to search for 18S rRNA instead of mtDNA, allow partial sequences, specify a blacklist of GIs and a minimum sequence length

`perl MConsensus.pl --out 18SNematoda --processors 20 --retmax 2000 --organism Nematoda --geneType rRNA --verbose --geneName "product\t18" --gene "18S+rRNA" --allowPartial --blackList 18SNematoda/blacklist.txt --minLen 100`

  Be careful when using the -allowPartial option, as this allows the consensus algorithm to keep all alignment positions, with no regard
  for how many sequences cover that position. It is a good idea to look at the consensusTable.txt file in excel and filter out
  any positions not covered by enough sequences.  