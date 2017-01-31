# MConsensus.pl

  Generate a consensus for a specified gene in a specified taxa

## A word of warning

  This script is a PITA. It requires a lot of babysitting of the output to make sure it's outputting something 
  meaningful and not a load of junk. If you're not comfortable looking through the perl code to figure things
  out I would recommend steering clear. If you do try to use it I would recommend at a minimum using the 
  --verbose option and watching the output. Also, look over your alignment output to be sure that it looks 
  like a decent alignment. If your alignment is bad, your consensus is going to be bad. There are a lot of 
  poorly annotated sequences on NCBI and they can mess up everything pretty badly in this script. In the 
  future I'd like to allow the script to accept a "black list" of accessions to omit from the analysis. 
  That way if you run it and find some bad sequences you can rerun it while ignoring them. But for now you 
  have to carefully watch the output. So, I would recommend not using this script unless you feel like wading
  into a bit of a mess potentially.
    
## perl MConsensus.pl [options]

  This script takes in an organism name (any scientific taxon) and a gene (or other genetic element) name.
  The script then downloads all the sequences for that query from the NCBI nt database. Those sequences are 
  then blasted against a local nt database to get any additional sequences (matching the input organism). 
  All the sequences are then combined and aligned using kalign. The consensus sequence of the alignment is
  then output as both a fasta and as a table of counts for each A/T/G/C at each position. 

## Dependencies

### This program is designed to work on Unix 
    
### blastx: blastn is used by this program and needs to be in your $PATH
    
  Be sure the $BLASTDB system variable is set to point to the directory containing the taxdb 
  database. This is done by export BLASTDB='absolute/path/to/your/blastdb/'

### kalign: kalign is used to align the sequences obtained and needs to be in your $PATH


## Output files

allSeqsAligned.fasta
- All sequences (both blast results and sequences retrieved from NCBI), aligned.

allSeqs.fasta
- All sequences (both blast results and sequences retrieved from NCBI).

blastFasta.fasta
- Fasta file of sequences from Blast results

blastResults.txt
- Results of Blasting the originalGis.fasta

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

--blastDb (nt)

  Location of the nt blast database.

--blastHitCount (10)

  Number of blast hits to retrieve per sequence retrieved from NCBI.

--processors (1)

  The number of processors to use.

--percIdentCutoff (80)

  This is the minimum percent similarity to be included in the blast results

--blastMinLen (100)

  This is the shortest sequence retrieved during the blast step that can be kept.

--minAF (20)

  This is the minimum allele frequency (in percentage) to include when generating the consensus sequence.
  A lower value will result in more ambiguous bases. 

--verbose

  This option prints out status updates while the program runs.

--help 

  Print out this "helpful" information. 
