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
my $gene = "mitochondr*+\"complete+genome\"";
my $email = "matthewvc1\@gmail.com";
my $geneNameToMatch = "gene\tcytb,gene\tcob";
my $geneType = "gene";
my $outDir = "";
my $kickOutDiffSeqs;
my $minIdentToKick = 20; # minimum percent diff from another seq to cast vote to kick
my $minPercentDiffToKick = 80; # minimum percent of kick votes to keep
my $maxInsertToKick = 10; # max number of inserts to vote to kick
my $minInsertPercentToKick = 70; # minimum percent of kick votes to keep
my $allowPartial;
my $blackList = "";
my $minLen = 0;
my $maxLen = "inf";
my $p = 1;
my $minAf = 20;
my $verbose;
my $help;

GetOptions ("retmax=i"                 => \$retmax,
            "organism=s"               => \$organism,
            "gene=s"                   => \$gene,
            "email=s"                  => \$email,
	    "geneNameToMatch=s"        => \$geneNameToMatch,
	    "geneType=s"               => \$geneType,
            "outDir=s"                 => \$outDir,
	    "kickOutDiffSeqs"          => \$kickOutDiffSeqs,
	    "minIdentToKick=i"         => \$minIdentToKick,
	    "minPercentDiffToKick=i"   => \$minPercentDiffToKick,
	    "maxInsertToKick=i"        => \$maxInsertToKick,
	    "minInsertPercentToKick=i" => \$minInsertPercentToKick,
	    "allowPartial"             => \$allowPartial,
	    "blacklist=s"              => \$blackList,
	    "minLen=i"                 => \$minLen,
	    "maxLen=i"                 => \$maxLen,
	    "processors=i"             => \$p,
	    "minAF=i"                  => \$minAf,
            "verbose"                  => \$verbose,
            "help"                     => \$help            
) or pod2usage(1) && exit;

pod2usage(2) && exit if ($help);


##############################
# Global variables
##############################
my $nucSearch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&retmode=json";
my @tempFiles;
my @giArray;
my %annotHash;
my $efetch = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=";
my $giRetrieveNum = 1;
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
my %blackHash;

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
close GITXT;

##############################
### Get fasta seqeunces for originalGis.txt and make tempdir/(originalGis.fasta)
### Include email in web query
### Report how long it took
my $seqResponse;
my @gisToGet = @giArray;
while(scalar(@gisToGet) > 200) {
    my @currentGis = splice(@gisToGet, 0, 200);
    my $command = "GET \"$efetch" . join(',', @currentGis) . "&email=". $email . ")\"";
    $seqResponse .= `$command`;
}
my $command = "GET \"$efetch" . join(',', @gisToGet) . "&email=". $email . ")\"";
$seqResponse .= `$command`;

open (GIFASTAS, ">", $outDir."originalGis.fasta") or die "Cannot create originalGis.txt, check permissions\n";
print GIFASTAS $seqResponse, "\n";
push @tempFiles, $outDir."originalGis.fasta";
if($verbose) {
    print STDERR "Fasta sequences downloaded and written to ", $outDir."originalGis.fasta\n\n";
}
close GIFASTAS;

##############################
### Get annotation for original sequences
@gisToGet = @giArray;
my $annotResponse;
my $annotCommand;
while(scalar(@gisToGet) > 200) {
    my @currentGis = splice(@gisToGet, 0, 200);
    $annotCommand = "GET \"eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=ft&retmode=text&id=" . join(',', @currentGis) . "&email=". $email . "\"";
    $annotResponse .= `$annotCommand`;
}
$annotCommand = "GET \"eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=ft&retmode=text&id=" . join(',', @gisToGet) . "&email=". $email . "\"";
$annotResponse .= `$annotCommand`;

open (GIANNOT, ">", $outDir."originalGisAnnot.txt") or die "Cannot create originalGisAnnot.txt, check permissions\n";
print GIANNOT $annotResponse, "\n";
push @tempFiles, $outDir."originalGisAnnot.fasta";
if($verbose) {
    print STDERR "Annotations downloaded and written to ", $outDir."originalGisAnnot.txt using $annotCommand\n\n";
}
close GIANNOT;
##############################
### Use annotation to modify fasta files to start at CYTB

$/ = "\n";
open (GIANNOT, "<", $outDir."originalGisAnnot.txt") or die "Cannot open originalGisAnnot.txt.\n";
my $lastline = "";
my $header;

my @matchGenes = split ",", $geneNameToMatch;
my $geneMatch = join("|", @matchGenes);

while(my $input = <GIANNOT>) {
    chomp $input;
    if($input =~ /^>/) {
        $header = $input;
	$header =~ s/.+?\|//;
	$header =~ s/\|.*//;
    }
    if($allowPartial) {
	$lastline =~ s/[\<\>]//g;
    }
    if($input =~ /$geneMatch/i && $lastline =~ /[0-9]+\t[0-9]+\t$geneType/i){
        my ($number, $number2, undef) = split "\t", $lastline;
	if($number !~ /[<>]/ && $number2 !~ /[<>]/) { 
	    if($number < $number2) {
		$annotHash{$header} = $number . "\t" . $number2;	# make hash of positions $hash{header} = pos;
	    } else {
		$annotHash{$header} = $number2 . "\t" . $number;
	    }
	}
	if($number eq "") {
	    print STDERR $header, "\t", $lastline, "\n";
	}
    }
$lastline = $input;
}
close GIANNOT;

##############################
### If provided, pull in blacklist of and create hash 
if($blackList ne "") {
    open (BLACKLIST, "$blackList") or die "Cannot open blacklist file\n";
    while(my $blackIn = <BLACKLIST>) {
	chomp $blackIn;
	$blackIn =~ s/[\s\t]//g;
	$blackHash{$blackIn} = 1;
    } 
close BLACKLIST;
}


##############################
### Use hash of positions to pull in each fasta entry and rearrange it, then print out fixed fasta
### Also fail to print any sequence found in blacklist
$/ = "\n>";
open (GIFASTA, "<", $outDir."originalGis.fasta") or die "Cannot open originalGis.fasta.\n";
open (GIFASTAFIXED, ">", $outDir."originalGisFixed.fasta") or die "Cannot open originalGisFixed.fasta.\n";
open (BADANNOT, ">", $outDir . "badAnnots.txt") or die "Cannot write to badAnnots.txt";
my $badCount = 0;
my %foundSpecies;
my $speciesCount = 0;
while(my $input = <GIFASTA>) {
    chomp $input;
    my ($header, @sequence) = split "\n", $input;
    $header =~ s/>//;
    my $shortHeader = $header;
    $shortHeader =~ s/\s.+//;
    my $species = $header;
    $species =~ s/.+?\s//;
    $species =~ s/\s/_/;
    $species =~ s/\s.+//;
    if(!defined($foundSpecies{$species})) {	
	if(defined($annotHash{$shortHeader}) && !defined($blackHash{$shortHeader})) {
	    my ($start, $end) = split "\t", $annotHash{$shortHeader};
	    my $seq = join("", @sequence);
	    $seq =~ s/[\t\s]//g;
	    my $fixedSeq = substr($seq, $start - 1, ($end - $start) + 1);
	    if(length($fixedSeq) > $minLen && length($fixedSeq) < $maxLen) { # get rid of any sequences too short or too long
		print GIFASTAFIXED ">", $header, "\n", $fixedSeq, "\n"; ########### Need to change this to keep the longest sequence
	    }
	    $speciesCount++;
	    $foundSpecies{$species} = 1;
	} else {
	    $badCount++;
	    print BADANNOT $header, "\n";
	}
	#$foundSpecies{$species} = 1;
	#$speciesCount++;
    } 
}
if($verbose && $badCount > 0) {
    print STDERR "A total of ", $badCount, " fasta annotations could not be parsed\nHeaders written to badAnnots.txt\nA total of $speciesCount unique species kept\n\n";
} elsif($verbose) {
    print STDERR "All of your entries could be parsed. Amazing!\nA total of $speciesCount unique species kept\n\n";
}
close GIFASTA;
close GIFASTAFIXED;
close BADANNOT;
$/ = "\n";

############################
### Align blastResults.fasta and originalGis.fasta (aligned.aln)

#system("cat " . $outDir . "blastFasta.fasta " . $outDir . "originalGisFixed.fasta > " . $outDir . "allSeqs.fasta");
my $alignCommand;


if($verbose) {
    print STDERR "Combining original fasta sequences and blast result sequences and aligning using mafft\n";
}
# Made mafft the default because of the ability to adjust the direction of the sequences
#$alignCommand = "mafft --adjustdirectionaccurately --globalpair --thread " . $p . " " . $outDir . "allSeqs.fasta > " . $outDir."allSeqsAligned.fasta";
$alignCommand = "mafft --quiet --maxiterate 1000 --localpair --adjustdirectionaccurately --thread " . $p . " " . $outDir . "originalGisFixed.fasta > " . $outDir."allSeqsAligned.fasta";
if($verbose) {
    $alignCommand =~ s/mafft --quiet/mafft/;
}

my $doit = `$alignCommand`;

if($verbose) {
    print STDERR "Alignment completed\n\n";
}

############################
### Read in alignment, count the number of differences between each sequence and kick out any over a given percent difference
### I want to kick out any sequence that is more than X% different from more than Y% of all sequences
### The intent is to remove any sequence that is completely off or has large insertions which can throw off the consensus

if($kickOutDiffSeqs) {
    my $kickedOutCount = 0;
    my %diffsSeqHash;

    if($verbose) {
	print STDERR "Removing sequences too different from all others\n";
	print STDERR "Check settings for this step: minIdentToKick:         $minIdentToKick\n";
	print STDERR "                              minPercentDiffToKick:   $minPercentDiffToKick\n";
	print STDERR "                              maxInsertToKick:        $maxInsertToKick\n";
	print STDERR "                              minPercentInsertToKick: $minInsertPercentToKick\n";
	print STDERR "Reading in alignment\n\n";
    }
    $/ = "\n>";
    open (ALIGN, "<", $outDir . "allSeqsAligned.fasta") or die "Cannot open input file\n";

    # Read in alignment by sequence
    while(my $input = <ALIGN>) {
	chomp $input;
	my($header, @sequences) = split "\n", $input;
	$header =~ s/^>//;                  # remove first ">"
	my $sequence = join("", @sequences);
	$sequence = uc($sequence);
	$diffsSeqHash{$header} = $sequence;
    }
    close ALIGN;
    my $command = "cp $outDir" . "allSeqsAligned.fasta $outDir" . "allSeqsAlignedOriginal.fasta"; # move the original alignment to save it
    my $doit = `$command`;

    open (DIFFCOUNTS, ">", $outDir . "alignedDiffCount.txt") or die "Cannot write to alignment difference count output file\n";

    if($verbose) {
	print STDERR "Filtering sequences\n";
    }

    my $kickedOutSeqsThisLoop = 1;
    my $loopCount = 1;
    my $kickedOutFromDiffs = 0;
    my $kickedOutFromInserts = 0;
    while($kickedOutSeqsThisLoop > 0) {
	my @keys = keys %diffsSeqHash;
	my @toRemove = (); # index of sequences to kick out
	for(my $i = 0; $i < scalar(@keys); $i++) {  # for each sequence, check if it should be kicked
	    my $tooDiffCount = 0;
	    my $highInsertCount = 0;
	    my $sequence = $diffsSeqHash{$keys[$i]};
	    for(my $j = 0; $j < scalar(@keys); $j++) {  # go through all the seqs
		my $compareSeq = $diffsSeqHash{$keys[$j]};
		my $diffCount = 0;
		my $insertCount = 0;
		
		# go through each base, count diffs and inserts
		for(my $k = 0; $k < (length($sequence) - 1); $k++) { 
		    my $seqBase = substr($sequence, $k, 1);
		    my $compareBase = substr($compareSeq, $k, 1);
		    if($seqBase ne $compareBase) {
			############################# make exception for insertions???
			$diffCount++;
			if($compareBase eq "-") { # if both are "-" it isn't an insert
			    $insertCount++;
			}
		    }
		}          
                # percent similarity
		if(100 * (1 - ($diffCount / length($sequence))) < $minIdentToKick) { # too different
		    $tooDiffCount++;
		}
		# percent inserts
		if(100 * ($insertCount / length($sequence)) > $maxInsertToKick) { # too many inserts
		    $highInsertCount++;
		}

		print DIFFCOUNTS $loopCount, "\t", $keys[$i], "\t", $keys[$j], "\t", 
		        100 * (1 - ($diffCount / length($sequence))), "\t",
            		100 * ($insertCount / length($sequence)), "\n";

	    }         #Percent good
	    if( 100 * (1 - ($tooDiffCount / scalar(@keys))) < $minPercentDiffToKick) {
		push @toRemove, $i;
		$kickedOutFromDiffs++;
	    } elsif(100 * ($highInsertCount / scalar(@keys)) > $minInsertPercentToKick) {
		push @toRemove, $i;
		$kickedOutFromInserts++;
	    }
	}
	for my $index (@toRemove) {
	    delete $diffsSeqHash{$keys[$index]};
	}
	$kickedOutSeqsThisLoop = scalar(@toRemove);
	$kickedOutCount += scalar(@toRemove);
	$loopCount++;
    }
    close DIFFCOUNTS;

    if($verbose) {
	print STDERR "Printing out new alignment and summary file\n";
    }

    # Print it out here
    open (NEWSEQS, ">", $outDir . "allSeqsFiltered.fasta") or die "Cannot write to new alignment file\n";
    for my $header (keys %diffsSeqHash) {
	my $seq = $diffsSeqHash{$header};
	$seq =~ s/-//g;
	print NEWSEQS ">", $header, "\n", $seq, "\n";
    }

    #open (NEWALIGN, ">", $outDir . "allSeqsAligned.fasta") or die "Cannot write to new alignment file\n";
    my $newAlignCommand = $alignCommand;
    $newAlignCommand =~ s/originalGisFixed.fasta/allSeqsFiltered.fasta/;
    my $doit2 = `$newAlignCommand`;

    if($verbose) {
        print STDERR "Realignment completed\n\n";
    }

  
    if($verbose) {
	print STDERR "Removed $kickedOutCount sequences due to excess differences. Diffs: $kickedOutFromDiffs. Inserts: $kickedOutFromInserts\n\n";
    }
}

$/ = "\n";





############################
# Parse aligned.aln with perl script to generate two things:
#       -consensus sequence based on minumum %RAF/MiAF (consensus.fasta)
#       -table file with one row per nt and one column per A/T/G/C with counts for
#         each nucleotide at each base and one column of consensus (seqCounts.txt)

if($verbose) {
    print STDERR "Generating consensus from ", $outDir . "allSeqsAligned.fasta\n";
}


$/  = "\n>"; # change input delimiter
open (ALIGN, "<", $outDir . "allSeqsAligned.fasta") or die "Cannot open input file\n";
while(my $input = <ALIGN>) {
    chomp $input;
    my($header, @sequences) = split "\n", $input;
    $header =~ s/^.+?\s//;                  # remove first ">"
    #$header =~ s/^[0-9]+//;             # remove gi info from original gis
    #$header =~ s/^gi.+\|\s//;           # get rid of gi info
    $header =~ s/\s/_/;                 # replace genus species with genus_species
    $header =~ s/\s.+//;                # remove the rest of the crap in species name (isolate, etc..)
    my $sequence = join("", @sequences);
    $sequence = uc($sequence);
    $sequence =~ s/[YRWSKMDVHBN]/-/g;   # replace ambiguous bases with a dash
    for(my $i = 0; $i < length($sequence); $i++) {
        my $base = substr($sequence, $i,1);
        $alignmentHash{$i}{$header}{$base}++;
    }
    if($verbose) {
        $seqCount++;
        print STDERR $seqCount,"\r";
    }
}

close ALIGN;
$/  = "\n"; # change input delimiter back

open (CONSENSUSFILE, ">", $outDir . "Consensus.fasta") or die "Cannot write to consensus output file\n";
open (TABLEFILE, ">", $outDir . "ConsensusTable.txt") or die "Cannot write to consensus table output file\n";

print TABLEFILE "AlignmentBase\tA\tT\tG\tC\tConsensus\tMissing\n";
print CONSENSUSFILE ">", $organism, "_", $gene, "_", $geneNameToMatch, "\n";

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
    my $countCutoff = ($minAf / 100) * $speciesCount; # minimum number of counts
    if($missingBase < $countCutoff || $allowPartial) {  # keep only those bases that have info at more than $minAf % bases 
	                                                # keep all bases if allow partial is enabled - otherwise the consensus can be super short 

	print TABLEFILE $baseNum, "\t";
	if($speciesCount > 0) {
	    #my $countCutoff = ($minAf / 100) * $speciesCount; # minimum number of counts
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
}

print CONSENSUSFILE "\n";
print STDERR "\n";
close CONSENSUSFILE;

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

=head NAME

    MConsensus.pl - generates a consensus for a specified gene in a specified taxa
    
=head SYNOPSIS
    
    perl MConsensus.pl [options]

=head OPTIONS

    This script takes in an organism name (any scientific taxon) and a gene (or other genetic element) name.
    The script then downloads all the sequences for that query from the NCBI nt database. Those sequences are 
    then blasted against a local nt database to get any additional sequences matching the input organism. 
    All sequences are then combined and aligned using mafft. The consensus sequence of the alignment is
    then output as both a fasta and as a table of counts for each A/T/G/C at each position. 

    The dependancies for this program are:
    
=over 4

=item B<mafft>

  mafft: mafft is used to align the sequences obtained and needs to be in your $PATH

=back

    
    The options used by this program are:

=over 4

=item B<--retmax> (1000)

    The maximum number of hits to retrieve from the NCBI database matching the organism and gene

=item B<--organism> ("salamanders")

    The organism to be queried. Use scientific names (e.g. canidae) or other names recognized by 
    NCBI database (e.g. salamanders). I would highly recommend testing your search here: 
        http://www.ncbi.nlm.nih.gov/nuccore/

=item B<--gene> ("mitochondr*+\"complete+genome\"")

    The gene sequence (or other genetic element) to be queried. If you use multiple words, put them 
    in quotes and put a plus between the words. I would highly recommend testing your search here: 
        I<http://www.ncbi.nlm.nih.gov/nuccore/>

=item B<--geneNameToMatch> ("gene\tcytb,gene\tcob")

    Fill this in later

=item B<--geneType> ("gene")

    Fill this in later

=item B<--email> (matthewvc1@gmail.com)

    Your email. This is used in the query to NCBI.

=item B<--outDir> (organism_date_time)

    The directory to write the files out to. By default, the directory is named for the organism and
    month/day/time the script is run.

=item B<--processors> (1)

    The number of processors to use.

=item B<--allowPartial> 

    This is a flag that tells the program to allow partial sequences to be used in the alignment and consensus.
    A partial sequence is defined as having either ">" or "<" in either of the gene position fields in the 
    annotation. 

=item B<--blackList>

    This is a simple list of accessions to exclude from the alignment and consensus. This is found after the
    ">" in the fasta file and up to (but not including) the first space. eg: KY119968.1. 

=item B<--minAF> (20)

    This is the minimum allele frequency (in percentage) to include when generating the consensus sequence.
    A lower value will result in more ambiguous bases. 

=item B<--minLen> (0)

  Minimum sequence length to be kept.

=item B<--maxLen> (INF)

  Maximum sequence length to be kept.

=item B<--verbose>

    This option prints out status updates while the program runs.

=item B<--help> 

    Print out this "helpful" information. 

=back

=cut
