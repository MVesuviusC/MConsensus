#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

##############################
# By Matt Cannon
# Date: May 7, 2016
# Last modified: June 10, 2016
# Title: MConsensus.pl
# Purpose: generate consensus sequence for a specified gene in a specified taxon
##############################

##############################
# Options
##############################


my $retmax = 1000;
my $organism = "salamanders";
my $gene = "mitochondrion+\"complete+genome\"";
my $email = "matthewvc1\@gmail.com";
my $outDir = "";
my $blastDb = "nt";
my $blastHitCount = 10;
my $p = 1;
my $clustalo;
my $mafft;
my $pidentCutoff = 80;
my $minAf = 20;
my $blMinLen = 100;
my $verbose;
my $help;

GetOptions ("retmax=i"          => \$retmax,
            "organism=s"        => \$organism,
            "gene=s"            => \$gene,
            "email=s"           => \$email,
            "outDir=s"          => \$outDir,
            "blastDb=s"         => \$blastDb,
            "blastHitCount=i"   => \$blastHitCount,
            "processors=i"      => \$p,
	    "clustalo"          => \$clustalo,
	    "mafft"             => \$mafft,
            "percIdentCutoff=i" => \$pidentCutoff,
            "blastMinLen=i"     => \$blMinLen,
            "minAF=i"           => \$minAf,
            "verbose"           => \$verbose,
            "help"              => \$help            
) or pod2usage(1) && exit;

pod2usage(2) && exit if ($help);


##############################
# Global variables
##############################
my $nucSearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmode=json";
my @tempFiles;
my @giArray;
my $efetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=";
my $blastOutfmt = "-outfmt \"6 sgi staxids length pident evalue bitscore sseq\"";
my $giRetrieveNum = 1;
my %blastResultsHash;
my %alignmentHash;
my %ambiguityHash = (
        A     => "A",
        T     => "T",
        C     => "C",
        G     => "G",
        TC    => "Y",
        AG    => "R",
        AT    => "W",
        GC    => "S",
        TG    => "K",
        AC    => "M",
        ATG   => "D",
        AGC   => "V",
        ATC   => "H",
        TGC   => "B", 
        ATGC  => "N",
        N     => "N");
my $seqCount = 0;
my @speciesList;

##############################
# Code
##############################

##############################
### Make an output directory (output[date/time]/)
### Make sure have write permission
if(-w "." eq "") { # check for write permissions
      print STDERR "No write permissions to this directory!\n\n";
      die;
}
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
if($outDir eq "") {
      $outDir = $organism . "M" . $mon . "_D" . $mday . "_H" . $hour . "_M" . $min . "_S" . $sec . "/";
} else {
      $outDir .= "/";
}
mkdir $outDir;
if($verbose) {
    print STDERR "Output directory: ", $outDir, "\n\n";
}

##############################
### Get gis from nucleotide database web query and put in tempdir/(originalGis.html)
my $esearch = "GET \"" . $nucSearch . "&retmax=" . $retmax . "&email=" . $email . "&term=(" . $organism . "[Organism]" . "+AND+" . $gene . ")\"";
if($verbose) {
    print STDERR "Retrieving list of gis from: ", $esearch, "\n";
}
my $response = `$esearch`;
open (GIHTML, ">", $outDir."originalGis.html") or die "Cannot open originalGis.html, check permissions\n";
print GIHTML $response, "\n";
push @tempFiles, $outDir."originalGis.html";
close GIHTML;
if($verbose) {
    print STDERR "Gis retrieved and printed to ", $outDir."originalGis.html\n\n";
}

##############################
### Parse esearch and return list of gis tempdir/(originalGis.txt)
### Report how many gis returned and how long it took
if($verbose) {
    print STDERR "Parsing gis html file\n";
}
open (GITXT, ">", $outDir."originalGis.txt") or die "Cannot create originalGis.txt, check permissions\n";

$response =~ s/[\n\s\"]//g;
# get rid of everything up to the id list
$response =~ s/.+idlist:\[//; 
$response =~ s/].+//;

@giArray = split ",", $response;
print GITXT join("\n", @giArray), "\n";
push @tempFiles, $outDir."originalGis.txt";
if($verbose) {
    print STDERR scalar(@giArray), " gis parsed and written to ", $outDir."originalGis.txt", "\n\n";
}
if(scalar(@giArray) == 0) {
    print STDERR "No sequences retrieved from NCBI. Please check your search terms at http://www.ncbi.nlm.nih.gov/nuccore/\n";
    die;
}

##############################
### Get fasta seqeunces for originalGis.txt and make tempdir/(originalGis.fasta)
### Include email in web query
### Report how long it took
my $command = "GET \"$efetch" . join(',', @giArray) . "&email=". $email . ")\"";
my $seqResponse = `$command`;
open (GIFASTAS, ">", $outDir."originalGis.fasta") or die "Cannot create originalGis.txt, check permissions\n";
print GIFASTAS $seqResponse, "\n";
push @tempFiles, $outDir."originalGis.fasta";
if($verbose) {
    print STDERR "Fasta sequences downloaded and written to ", $outDir."originalGis.fasta\n\n";
}


##############################
### Generate list of gis for $organism to use in blast search
### Report how long it took

my $giSearch = "GET \"" . $nucSearch . "&retmax=100000" . "&email=" . $email . "&term=(" . $organism . "[Organism]" . ")\"";
if($verbose) {
    print STDERR "Retrieving list of gis matching ", $organism, " from: ", $giSearch, "\n";
    print STDERR "Retrieving batch #", $giRetrieveNum, " from ncbi.\n";
}
my $giResponse = `$giSearch`;

$giResponse =~ s/[\n\s\"]//g;
# get rid of everything up to the id list
$giResponse =~ s/.+idlist:\[//; 
$giResponse =~ s/].+//;
my @orgGiArray = split ",", $giResponse;

while(scalar(@orgGiArray % 100000 == 0)) {
    if($verbose) {
        print STDERR "Retrieving batch #", $giRetrieveNum + 1, " from ncbi.\n";
    }
# get more gis, starting after the first 100,000
    $giSearch = "GET \"" . $nucSearch . "&retmax=100000&retstart=" . $giRetrieveNum * 100000 . "&email=" . $email . "&term=(" . $organism . "[Organism]" . ")\"";
$giResponse = `$giSearch`;
    $giResponse =~ s/[\n\s\"]//g;
    # get rid of everything up to the id list
    $giResponse =~ s/.+idlist:\[//; 
    $giResponse =~ s/].+//;
    my @tempGiArray = split ",", $giResponse;
    push @orgGiArray, @tempGiArray; 
    $giRetrieveNum++;
}

open (ORGGITXT, ">", $outDir."organismGis.txt") or die "Cannot create organismGis.txt, check permissions\n";
print ORGGITXT join("\n", @orgGiArray), "\n";
push @tempFiles, $outDir."organismGis.txt";
if($verbose) {
    print STDERR scalar(@orgGiArray), " gis parsed and written to ", $outDir."organismGis.txt", "\n\n";
}


############################
### Blast originalGis.fasta against taxaDb and write out tempdir/blastResults.txt
### Keep $blastHitCount hits per query
### blastResults.txt should use outfmt 6 with columns:
###       - gi(sgi), taxaInfo(staxids), alignment length(length), 
###         percent identical matches(pident), evalue(evalue), 
###         bit score(bitscore), scientific name(sscinames), aligned portion of subject(sseq)
### Report how long it took 

if(scalar(@giArray) > 0) {
    if($verbose) {
     print STDERR "Blasting originalGis.fasta against ", $blastDb, "\n";
    }
    system("blastn", 
        "-gilist", $outDir."organismGis.txt", 
        "-query", $outDir."originalGis.fasta", 
        "-db", $blastDb,  
        "-out", $outDir."blastResults.txt",
        "-outfmt", "6 sgi staxids length pident evalue bitscore sscinames sseq",
        "-num_threads", $p, 
        "-num_alignments", $blastHitCount );
    if($verbose) {
     print STDERR "Blast finished.Output written to ",$outDir."blastResults.txt", "\n\n";
    }
} else {
    print STDERR "No sequences matching search terms retrieved.\nCheck search terms at https://www.ncbi.nlm.nih.gov/nuccore/\n";
    die;
}


############################
### Parse blastResults.txt 
###       - Keep one longest sequence per species
###       - Cutoffs for minimum length and % identity
### Write fasta file (blastResults.fasta) with sequence and header: >gi_species

# What about two blast hits that hit different parts of the query, but are from the same species?

if($verbose) {
    print STDERR "Parsing blast results\n";
}
open (BLASTRESULTS, "<", $outDir."blastResults.txt") or die "Cannot open blastResults.txt\n";
while(my $blastInput = <BLASTRESULTS>){
    chomp $blastInput;
    my ($sgi, $staxids, $length, $pident, $evalue, $bitscore, $sscinames, $sseq) = split "\t", $blastInput;
    $sseq =~ s/-//g; # get rid of insertions in the alignment
    if(exists($blastResultsHash{$sscinames})){
        if(length($blastResultsHash{$sscinames . " " . $sgi}) < length($sseq) && $pident >= $pidentCutoff) {
            $blastResultsHash{$sscinames . " " . $sgi} = $sseq;
        }
    } elsif($pident >= $pidentCutoff && $length >= $blMinLen) {
        $blastResultsHash{$sscinames . " " . $sgi} = $sseq;
    }
}

open (BLASTFASTA, ">", $outDir."blastFasta.fasta") or die "Cannot write to blastFasta.fasta\n";
for my $header (keys %blastResultsHash) {
    print BLASTFASTA ">", $header, "\n", $blastResultsHash{$header}, "\n";
}
if($verbose) {
    print STDERR "Blast results parsed. Output written to ", $outDir."blastFasta.fasta\n\n";
}

############################
### Align blastResults.fasta and originalGis.fasta (aligned.aln)

system("cat " . $outDir . "blastFasta.fasta " . $outDir . "originalGis.fasta > " . $outDir . "allSeqs.fasta");
my $alignCommand;

if($clustalo) {
    # Clustalo was slower
    if($verbose) {
	print STDERR "Combining original fasta sequences and blast result sequences and aligning using clustalo\n";
    }
    $alignCommand = "clustalo" . " --threads " . $p . " -i " . $outDir . "allSeqs.fasta" . " --out " . $outDir . "allSeqsAligned.fasta";
} elsif($mafft) {
    if($verbose) {
	print STDERR "Combining original fasta sequences and blast result sequences and aligning using mafft\n";
    }
    $alignCommand = "mafft --adjustdirectionaccurately --globalpair --thread " . $p . " " . $outDir . "allSeqs.fasta > " . $outDir."allSeqsAligned.fasta";

} else {
    if($verbose) {
	print STDERR "Combining original fasta sequences and blast result sequences and aligning using kalign\n";
    }
    $alignCommand = "kalign -i " . $outDir . "allSeqs.fasta" . " -o " . $outDir . "allSeqsAligned.fasta";
}
my $doit = `$alignCommand`;

if($verbose) {
    print STDERR "Alignment completed\n\n";
}

############################
# Parse aligned.aln with perl script to generate two things:
#       -consensus sequence based on minumum %RAF/MiAF (consensus.fasta)
#       -table file with one row per nt and one column per A/T/G/C with counts for
#         each nucleotide at each base and one column of consensus (seqCounts.txt)

if($verbose) {
    print STDERR "Generating consensus from ", $outDir . "allSeqsAligned.fasta\n";
}

#system("perl scripts/summarizeAlignment.pl --input " . $outDir."allSeqsAligned.fasta");

$/  = "\n>"; # change input delimiter
open (ALIGN, "<", $outDir . "allSeqsAligned.fasta") or die "Cannot open input file\n";
while(my $input = <ALIGN>) {
    chomp $input;
    my($header, @sequences) = split "\n", $input;
    $header =~ s/^>//;                  # remove first ">"
    $header =~ s/^[0-9]+//;             # remove gi info from original gis
    $header =~ s/^gi.+\|\s//;           # get rid of gi info
    $header =~ s/\s/_/;                 # replace genus species with genus_species
    $header =~ s/\s.+//;                # remove the rest of the crap in species name (isolate, etc..)
    my $sequence = join("", @sequences);
    $sequence =~ s/[YRWSKMDVHBN]/-/g;   # replace ambiguous bases with a space
    for(my $i = 0; $i < length($sequence); $i++) {
        my $base = substr($sequence, $i,1);
        $alignmentHash{$i}{$header}{$base}++;
    }
    if($verbose) {
        $seqCount++;
        print STDERR $seqCount,"\r";
    }
}
$/  = "\n"; # change input delimiter back

open (CONSENSUSFILE, ">", $outDir . "Consensus.fasta") or die "Cannot write to consensus output file\n";
open (TABLEFILE, ">", $outDir . "ConsensusTable.txt") or die "Cannot write to consensus table output file\n";

print TABLEFILE "Base\tA\tT\tG\tC\tConsensus\tMissing\n";
print CONSENSUSFILE ">$organism\n";

for my $baseNum ( sort {$a <=> $b} keys %alignmentHash) {
    my %currentBase;
    my $speciesCount = 0;
    my @baseArray = ("A", "T", "G", "C");
    my $missingBase = 0;
    my $consensusBase = "";
    for my $species (keys %{ $alignmentHash{$baseNum} }) {      # Generate a consensus sequence per species
        if($baseNum == 1) {
            push @speciesList, $species;
        }
        my $topBase = "-";
        my $topBaseCount = 0;
        for my $base (keys %{ $alignmentHash{$baseNum}{$species} }) {   # Loop through the sequence
            if($base ne "-" && $alignmentHash{$baseNum}{$species}{$base} > $topBaseCount) {
                $topBase = $base;                               # keep the most common base for each position
                $topBaseCount = $alignmentHash{$baseNum}{$species}{$base};
            } elsif ($base ne "-" && $alignmentHash{$baseNum}{$species}{$base} == $topBaseCount && rand() >= 0.5 ) { 
                # 50% chance to replace current base if two+ bases found
                $topBase = $base;
            }
        }
        if($topBase ne "-") {
            $currentBase{$topBase}++;
            $speciesCount++;
        } else {
            $missingBase++;
        }
    }
    print TABLEFILE $baseNum, "\t";
    if($speciesCount > 0) {
        my $countCutoff = ($minAf / 100) * $speciesCount; # minimum number of counts
        for my $topBase (@baseArray) {
            if(exists($currentBase{$topBase})) {
                print TABLEFILE $currentBase{$topBase}, "\t";
                if($currentBase{$topBase} >= $countCutoff) {
                    $consensusBase .= $topBase;
                }
            } else {
                print TABLEFILE "0\t";
            }
        }
        if ( $consensusBase eq "" ) {   # If all bases are below the min %cutoff
            $consensusBase = "N";
        }
        print TABLEFILE $ambiguityHash{$consensusBase}, "\t";
    } else {
        print TABLEFILE "0\t0\t0\t0\t";
        $consensusBase = "N";
    }
    print TABLEFILE $missingBase, "\n";
    print CONSENSUSFILE $ambiguityHash{$consensusBase};
}

print CONSENSUSFILE "\n";
print STDERR "\n";

if($verbose) {
    @speciesList = sort(@speciesList);
    print STDERR "\n", join("\t", @speciesList), "\n";
}

if($verbose) {
    print STDERR "Consensus generated. Output written to ", $outDir."Consensus.fasta and ", $outDir."ConsensusTable.txt\n";
}


##############################
# Generate a .Rmd file (plot.Rmd) so I can generate a plot from seqCounts.txt with 
# bases along the x axis and a stacked bar plot for each base showing %of each nucleotide
# Also generate another plot of stacked bar plots with RAF and MiAF







##############################
# POD
##############################

    
=head SYNOPSIS

   
    MConsensus.pl - generates a consensus for a specified gene in a specified taxa
    
    perl MConsensus.pl [options]

=head OPTIONS

    This script takes in an organism name (any scientific taxon) and a gene (or other genetic element) name.
    The script then downloads all the sequences for that query from the NCBI nt database. Those sequences are 
    then blasted against a local nt database to get any additional sequences matching the input organism. 
    All sequences are then combined and aligned using kalign. The consensus sequence of the alignment is
    then output as both a fasta and as a table of counts for each A/T/G/C at each position. 

    The dependancies for this program are:
    
=over 4

=item B<blastx>
    
  blastx: blastn is used by this program and needs to be in your $PATH
    Be sure the $BLASTDB system variable is set to point to the directory containing the taxdb 
    database. This is done by export BLASTDB='absolute/path/to/your/blastdb/'

=item B<kalign>

  kalign: kalign is used to align the sequences obtained and needs to be in your $PATH

=back

    I need to write an explanation of what the output files are.


    The options used by this program are:

=over 4

=item B<--retmax> (1000)

    The maximum number of hits to retrieve from the NCBI database matching the organism and gene

=item B<--organism> ("salamanders")

    The organism to be queried. Use scientific names (e.g. canidae) or other names recognized by 
    NCBI database (e.g. salamanders). I would highly recommend testing your search here: 
        I<http://www.ncbi.nlm.nih.gov/nuccore/>

=item B<--gene> ("mitochondria+\"complete+genome\"")

    The gene sequence (or other genetic element) to be queried. If you use multiple words, put them 
    in quotes and put a plus between the words. I would highly recommend testing your search here: 
        I<http://www.ncbi.nlm.nih.gov/nuccore/>

=item B<--email> (matthewvc1@gmail.com)

    Your email. This is used in the query to NCBI.

=item B<--outDir> (organism_date_time)

    The directory to write the files out to. By default, the directory is named for the organism and
    month/day/time the script is run.

=item B<--blastDb> (nt)

    Location of the nt blast database.

=item B<--blastHitCount> (10)

    Number of blast hits to retrieve per sequence retrieved from NCBI.

=item B<--processors> (1)

    The number of processors to use.

=item B<--percIdentCutoff> (80)

    This is the minimum percent similarity to be included in the blast results

=item B<--blastMinLen> (100)

    This is the shortest sequence retrieved during the blast step that can be kept.

=item B<--minAF> (20)

    This is the minimum allele frequency (in percentage) to include when generating the consensus sequence.
    A lower value will result in more ambiguous bases. 

=item B<--verbose>

    This option prints out status updates while the program runs.

=item B<--help> 

    Print out this "helpful" information. 

=back

=cut
