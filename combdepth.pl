#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_compat permute no_getopt_compat);
#I am not activating bundling since it is not compatible with array inputs, otherwise I would use gnu_getopt
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::File;
use List::Util qw(min sum);
#use Data::Dumper;

## Conf
############
my $MAX_IT=200; ##Increase this if working with contings/scaffolds instead of chromosomes

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


# Initializing filehandles
my @files;
for (my $nfile=0; $nfile<scalar @inputfiles; ++$nfile)
{
	my $filehandle;

	! -s $inputfiles[$nfile] and die "The input file $inputfiles[$nfile] is empty!\n";

	if($inputfiles[$nfile] =~ /.gz$/)
	{
		$filehandle =  new IO::Uncompress::Gunzip $inputfiles[$nfile] or die "The file $inputfiles[$nfile] was detected as compressed with gzip, but it cannot be opened: $GunzipError\n";
	}
	else
	{
		$filehandle = IO::File->new($inputfiles[$nfile], "r") or die "The file $inputfiles[$nfile] was detected as a regular file, but it cannot be opened. Add the .gz extension if it has been compressed with gzip\n";
	}
	$files[$nfile]=$filehandle;
}

# Initializing filters
my @filters;
my $i=0;
for (my $val=$min; $val<=$max; $val=$val+$by)
{
	$filters[$i]=$val;
	++$i;
}
$filters[$i-1] != $max and push(@filters,$max);
#print("DEBUG:".join(",",@filters)."\n");

# Constants
my $nfiles= scalar @inputfiles;
my $nfilters= scalar @filters;

## The two categories below would be better implemented as the different attributes of an object for each file file.chr, file.pos, file.depth, file.filehandler... I do not want to spend time with the perl's weird OOP

# Traking variables.
my @cchr; #Strings, name of the current chromosome in the file i
my @cpos; #Integer, pos in the current chromosome in the file i
my @cdepth; #Integer, depth in the current pos in the current chrom in the file i
my @cnewchr; #Bool, file i is situated in a new chromosome

# Slice variables. Arrays with the index of files that have the current chr or pos. Used to function as masks (from a functional perspective)
#my @chrslice;
#my @posslice;
#my @validfilesslice;
#my @nextvalidfilesslice;

# Mask variables. Arrays with booleans in the index of files that have the current chr or pos. Array-based masks
my @chrmask;
my @posmask;
my @posvalidfilesmask;
my @validfilesmask;
my @nextvalidfilesmask;

# Control variables
my $stop=0;
$i=0;
my $j=0;
my $chr=1;
my $pos=0;
my $ifilt=0;
my $isample=0;
my $nvalidfiles=0;
my $npassthisfilter=0;

# Content variables
my $line;
my @results;

##Initialization rutines

#Results init
for ($i=0; $i <= $nfiles; ++$i)
{
	for ($j=0; $j< $nfilters; ++$j)
	{
		$results[$i][$j]=0;
	}
}

#Loading first data point
for (my $nfile=0; $nfile<$nfiles; ++$nfile)
{
	$line=$files[$nfile]->getline();
	chomp($line);
	($cchr[$nfile],$cpos[$nfile],$cdepth[$nfile])=split("\t",$line);
	#print("DEBUG: chr=$cchr[$nfile], pos=$cpos[$nfile], depth=$cdepth[$nfile]\n");
}

$i=0;

##Main loop that gets new chromosomes until the end of all files or an iterator limit (to avoid infinite loops with a while)
while (!($stop == 1 || $i>=$MAX_IT))
{
	$pos=min(@cpos);

	#@chrslice=makeslice(\@cchr,$chr);
	#scalar @chrslice or die "Different chromosomes not implemented yet\n";
	#@chrmask=makemask(\@cchr,\@chrslice);

	@chrmask=makemask(\@cchr,$chr);
	sum @chrmask or die "Different chromosomes not implemented yet\n";

	#@posslice=makeslice(\@cpos,$pos);
	#scalar @posslice or die "None of the files contains this position, still implementing this\n";
	#@posmask=makemask(\@cpos,\@posslice);
	
	@posmask=makemask(\@cpos,$pos);
	sum @posmask or die "None of the files contains this position, still implementing this\n";
	
	@validfilesmask=arrayfunc(\@posmask, \&vand, \@chrmask);
	@posvalidfilesmask=@validfilesmask;
	@nextvalidfilesmask=@validfilesmask;
	$nvalidfiles=sum @validfilesmask;

	#@validfilesslice=makeslice(\@validfilesmask,1);	
	#$nvalidfiles=scalar @validfilesslice;
	#@nextvalidfilesslice=();

	for ($ifilt=0; $ifilt< scalar @filters; ++$ifilt) #desired depth
	{
		$npassthisfilter=0;
		for($isample=0; $isample<$nvalidfiles; ++$isample)
		{
			if($validfilesmask[$isample]==1)
			{
				if($cdepth[$isample]>=$filters[$ifilt])
				{
		#foreach $isample (@validfilesslice) #numberid of valid samples in this position ##TODO: I may need to evaluate if I should use masks with ifs or slices with the push to generate the nextslice. Which one will be quicker?
		#{
		#if($cdepth[$isample]>=$filters[$ifilt])
		#{
		#push(@nextvalidfiles,$isample); #this sample has enough coverage at this level and should be checked at the next
				
					++$npassthisfilter;
				}
				else
				{
					$nextvalidfilesmask[$isample]=0; #this sample does not have enough coverage at this level and therefore it does not need to be checked again for this position
					--$nvalidfiles;
				}
			}
		}
		++$results[$npassthisfilter][$ifilt];
		
		#@validfilesslice=@nextvalidfilesslice; #we discard the files that did not meet this filtering level
		#$nvalidfiles=scalar @nextvalidfilesslice;
		#@nextvalidfilesslice=();
		@validfilesmask=@nextvalidfilesmask;
	}
		
		
	#updatepos();	
	++$i;
	last;
}

# Print results
# #############

my @header;

push(@header,">=$_") foreach @filters;

print(join(",","Sharedby",@header),"\n");
for ($i=1; $i <= $nfiles; ++$i)
{
	for ($j=0; $j< scalar @filters; ++$j)
	{
		print(join(",",$i,@{$results[$i]}),"\n");
	}
}

## Functions
############

sub makemask
{
	my $n=scalar @{$_[0]};
	my @mask=(0) x $n;
	for (my $i=0; $i<$n; ++$i)
	{
		if(${$_[0]}[$i]==$_[1])
		{
			$mask[$i]=1;
		}
	}
	return @mask;
}

sub arrayfunc
{
	my @out;	
	
	for (my $j=0;$j<scalar @{$_[0]};++$j)
	{
		$out[$j]=$_[1]->($_[0]->[$j],$_[2]->[$j]);
	}
	
	return @out;
}

sub vand
{
	if (scalar @_ !=2)
	{
		return undef;
	}
	else
	{
		return $_[0] & $_[1];
	}	
}

sub vor
{
	if (scalar @_ !=2)
	{
		return undef;
	}
	else
	{
		return $_[0] | $_[1];
	}	
}
#sub makeslice
#{
#	my @out;
#	for (my $i=0; $i<$nfiles; ++$i)
#	{
#		if(${$_[0]}[$i]==$_[1])
#		{
#			push(@out,$i);
#		}
#	}
#	return @out;
#}
#
#sub makemask
#{
#	my @mask=(0) x scalar @{$_[0]};
#	@mask[@{$_[1]}]=(1) x scalar @{$_[1]};
#	return @mask;
#}
