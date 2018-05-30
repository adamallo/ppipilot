#!/opt/local/bin/perl
###!/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_compat permute no_getopt_compat);
#I am not activating bundling since it is not compatible with array inputs, otherwise I would use gnu_getopt

## Getopt I/O
#############

my @inputfiles;
my $output_prefix="";
my $min=10;
my $max=50;
my $by=10;
my $help;

my $usage="Usage: $0 -i inputsample1.tsv [inputsample2.tsv ... inputtsamplen.tsv] -l min -r max -b by -o output_prefix\n\nThis script takes as input the output of a number of \"samtools depth\" runs and calculates the number of common nucleotide positions at n= (max-min)/by depths in a range [min,max]\n";

(! GetOptions(
    'input|i=s{1,}' => \@inputfiles,
	'output|o=s' => \$output_prefix,
    'min|l=i' => \$min,
    'max|r=i' => \$max,
    'by|b=i' => \$by,
    'help|h' => \$help,
                )) or (($output_prefix eq "") || $help) and die $usage;


