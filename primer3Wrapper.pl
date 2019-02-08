#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

##############################
# By Matt Cannon
# Date:
# Last modified:
# Title: primer3Wrapper.pl
# Purpose: generate primer sequences using primer3
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $inputFasta;
my $minLen = 200;
my $maxLen = 450;
my $optSize = 350;
my $minTm = 57;
my $optTm = 60;
my $maxTm = 63;
my $maxAmbBases = 2;
my $misPrimingFasta;
my $primerName;
my $messy;
my $tile;
my $tileSize = 200;
my $noHeader;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
            "inputFasta=s"	=> \$inputFasta,
            "minLen=i"		=> \$minLen,
            "maxlen=i"		=> \$maxLen,
	    "optSize=i"         => \$optSize,
	    "minTm=i"           => \$minTm,
	    "optTm=i"           => \$optTm,
	    "maxTm=i"           => \$maxTm,
            "maxAmbBases=i"	=> \$maxAmbBases,
            "misPrimingFasta=s"	=> \$misPrimingFasta,
	    "primerName=s"      => \$primerName,
	    "leaveAMess"        => \$messy,
	    "tile"              => \$tile,
	    "tileSize=i"        => \$tileSize,
	    "noHeader"          => \$noHeader
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################
my $seq;
my $primer3InputText;
my $randomNum = int(rand(10000));
my @leftPrimerLocs;
my @rightPrimerLocs;
my $minSets = 1;
my @leftPrimers;
my @leftPrimersStart;
my @rightPrimers;
my @rightPrimersStart;

##############################
# Code
##############################


##############################
### Get input fasta
$/ = undef;
open INPUTFASTA, $inputFasta or die "cannot open input fasta\n";
my $input = <INPUTFASTA>;
# Check if the input is a fasta
if ($input !~ /^>/) {
    die "fasta input not in fasta format!\n";
}

my ($header, @sequence) = split "\n", $input;
$seq = join("", @sequence);
$seq =~ s/[\t\s]//g;

$/ = "\n"; 

# if $tile flag is defined, return 
if($tile) {
    $minSets = int(length($seq) / $tileSize);
}

if(!defined($primerName)) {
    $header =~ s/^>//;
    $primerName = $header;
}

# Check Tm settings to make sure they're not dumb
if($maxTm <= $minTm || $optTm < $minTm || $optTm > $maxTm) {
    print "Check your annealing Tm settings!\n";
    print "min Tm: $minTm, opt Tm: $optTm, max Tm: $maxTm\n";
    die;
}

for(my $i = 1; $i <= $minSets; $i++) {
    ##############################
    ### Make primer3 input file
    $primer3InputText = 
	"SEQUENCE_ID=" . $primerName . "\n" .
	"PRIMER_EXPLAIN_FLAG=1" . "\n" . 
	"SEQUENCE_TEMPLATE=" . $seq . "\n" . 
	"PRIMER_PRODUCT_SIZE_RANGE=" . $minLen . "-" . $maxLen . "\n" .
	"PRIMER_PRODUCT_OPT_SIZE=" . $optSize . "\n" .
	"PRIMER_MIN_TM=" . $minTm . "\n" .
	"PRIMER_OPT_TM=" . $optTm . "\n" .
	"PRIMER_MAX_TM=" . $maxTm . "\n" . 
	"PRIMER_LIBERAL_BASE=1\n";
    
    if(defined($misPrimingFasta)) {
	$primer3InputText .= "PRIMER_MISPRIMING_LIBRARY=" . $misPrimingFasta . "\n";
    }
    
    if($maxAmbBases > 0) {
	$primer3InputText .= "PRIMER_MAX_NS_ACCEPTED=" . $maxAmbBases . "\n";
    }
    
    if($tile) {
	$primer3InputText .= "SEQUENCE_TARGET=" . $i * $tileSize . ",10\n";
    }

    $primer3InputText .= "=\n";
    
    open PRIMER3FILE, ">", "primer3Input_" . $randomNum . ".txt";
    print PRIMER3FILE $primer3InputText;


    ##############################
    ### run primer3
    my $command = "primer3_core primer3Input_" . $randomNum . ".txt";
    my $primer3Output = `$command`;
    if(!$messy) {
	system("rm primer3Input_" . $randomNum . ".txt");
    }
    #print $primer3Output, "\n";

    ##############################
    ### Parse primer3 output

    my @primer3 = split "\n", $primer3Output;
    
    push @leftPrimers, grep(/PRIMER_LEFT_[0-9]+_SEQUENCE/, @primer3);
    s/PRIMER_LEFT_[0-9]+_SEQUENCE=// for @leftPrimers;
    
    push @leftPrimersStart, grep(/PRIMER_LEFT_[0-9]+=/, @primer3);
    s/PRIMER_LEFT_[0-9]+=// for @leftPrimersStart;
    s/,.+// for @leftPrimersStart;
    
    push @rightPrimers, grep(/PRIMER_RIGHT_[0-9]+_SEQUENCE/, @primer3);
    s/PRIMER_RIGHT_[0-9]+_SEQUENCE=// for @rightPrimers;
   
    push @rightPrimersStart, grep(/PRIMER_RIGHT_[0-9]+=/, @primer3);
    s/PRIMER_RIGHT_[0-9]+=// for @rightPrimersStart;
    s/,.+// for @rightPrimersStart;
}

if(!$noHeader) {
    print "FprimerName\tprimerF\tRprimerName\tprimerR\n";
}

if(scalar(@leftPrimers) == 0) {
    die "No primers found for $primerName\n";
}

for(my $i = 0; $i < scalar(@leftPrimers); $i++) {
    my $good = 1;
    my $j = $i - 1;
    if($i == 0) {
	print 
	    $primerName . "_" . $leftPrimersStart[$i] . "_F", 
	    "\t", $leftPrimers[$i], 
#	    "\t", $primerName . "_" . $rightPrimersStart[$i] . "_R",
	    "\t", $rightPrimersStart[$i] . "_R",
	    "\t", $rightPrimers[$i], "\n";
    } else {
	while($good == 1 && $j >= 0) {
	    # check if primer set is within 20bp of an already known pair
	    if(
	       abs($leftPrimersStart[$i] - $leftPrimersStart[$j]) < 20 &&
	       abs($rightPrimersStart[$i] - $rightPrimersStart[$j]) < 20 
	       ) {
		$good = 0;
	    }

	    $j--;

	}
	
	if($good == 1) {
	    print 
		$primerName . "_" . $leftPrimersStart[$i] . "_F", 
		"\t", $leftPrimers[$i], 
#		"\t", $primerName . "_" . $rightPrimersStart[$i] . "_R",
		"\t", $rightPrimersStart[$i] . "_R",
		"\t", $rightPrimers[$i], "\n";
	}
    }
}



##############################
# POD
##############################


=head SYNOPSIS

Summary:

    primer3Wrapper.pl - run primer3 on a fasta with options

Usage:

    perl primer3Wrapper.pl [options]


=head OPTIONS

Options:

    --verbose
    --help

=cut
