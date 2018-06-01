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

#Packages
#########

{
	package slurplikeline;
	use File::Slurp;
	use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

	sub new 
	{
		my ($class,$filename,$mode) = @_;
		my $n=0;
		my $maxi=0;
		my @lines;

		if (defined $mode && lc $mode eq 'gz')
		{
			gunzip $filename => \@lines;
		}
		else
		{
			#@lines = read_file($filename); #Not working
			my $filehandler;
			open($filehandler,$filename);
			@lines=<$filehandler>;
			close($filehandler);
		}

		$maxi=scalar @lines;
		
		my $self = bless { 
			filename => $filename,
			n => $n,
			lines => \@lines,
			maxi => $maxi
		},$class;

		return $self;
	}
	
	sub getline
	{
		my $self=shift;
		
		if ($self->{n} == $self->{maxi}){return undef};

		++$self->{n};
		return $self->{lines}->[$self->{n}-1];
	}
}



## Conf
############
my $MAX_IT=400000;
#my $MAX_IT=4000000000;

## Getopt I/O
#############

my @inputfiles;
my $output_prefix="";
my $min=10;
my $max=50;
my $by=10;
my $ref="";
my $dict="";
my $listchrs="";
my $help;
my $slurp=0;
my $yesall=0;

my $usage="Usage: $0 -i inputsample1.tsv [inputsample2.tsv ... inputtsamplen.tsv] --ref reference_genome.fa | --dict reference_genome.dict -o output_prefix [options]\nOptions: \n\t-l min \n\t-r max \n\t-b by \n\t-s read all files in memory (slurp). Much quicker but RAM intensive\n\t-y yes to all to use in combination with -s in non-interactive mode\n\t--list list of chromosmes sorted in the same order as the input tsv files (to use instead of ref or dict)\nThis script takes as input the output of a number of \"samtools depth\" runs and calculates the number of common nucleotide positions at n= (max-min)/by depths in a range [min,max]. \n";

(! GetOptions(
	'input|i=s{1,}' => \@inputfiles,
	'output|o=s' => \$output_prefix,
	'min|l=i' => \$min,
	'max|r=i' => \$max,
	'by|b=i' => \$by,
	'ref=s' => \$ref,
	'dict=s' => \$dict,
	'list=s' => \$listchrs,
	'slurp|s' => \$slurp,
	'yesal|y' => \$yesall,
	'help|h' => \$help,
									)) or (($output_prefix eq "") || $help) and die $usage;


# Measuring filesize and asking the user

if ($slurp == 1 && $yesall ==0)
{
	my $filesize;
	for (my $nfile=0; $nfile<scalar @inputfiles; ++$nfile)
	{
		$filesize+=-s $inputfiles[$nfile];
	}
	$slurp=prompt_yn("Are you sure that you want to slurp the files? This will use at least " . sprintf("%.2f",$filesize/1073741824) . "GB of RAM in your computer (it may be up to ".sprintf("%.2f",$filesize*6/1073741824)." GB if all files have been gzipped\n","n");
}

# Initializing filehandles
my @files;
for (my $nfile=0; $nfile<scalar @inputfiles; ++$nfile)
{
	my $filehandle;

	! -s $inputfiles[$nfile] and die "The input file $inputfiles[$nfile] does not exist or it is empty!\n$usage";

	if($inputfiles[$nfile] =~ /.gz$/)
	{
		if ($slurp == 1)
		{
			$filehandle = slurplikeline->new($inputfiles[$nfile],"gz") or die "Test";
		}
		else
		{
			$filehandle =  new IO::Uncompress::Gunzip $inputfiles[$nfile] or die "The file $inputfiles[$nfile] was detected as compressed with gzip, but it cannot be opened: $GunzipError\n";
		}
	}
	else
	{
		if ($slurp ==1)
		{
			$filehandle = slurplikeline->new($inputfiles[$nfile]) or die "Test";
		}
		else
		{
			$filehandle = IO::File->new($inputfiles[$nfile], 'r' ) or die "The file $inputfiles[$nfile] was detected as a regular file, but it cannot be opened. Add the .gz extension if it has been compressed with gzip\n";
			
		}
	}
	$files[$nfile]=$filehandle;
}

# Getting the list of chromosomes
my @chrs;
my %chrs;

if ($ref ne "")
{
	$dict=$ref;
	$dict=~s/\.fa/\.dict/;
}
if($dict ne "")
{
	open(my $filedict,$dict);
	@chrs=grep(/\@SQ/,<$filedict>);
	foreach my $ichr (@chrs)
	{
		$ichr=~s/\@SQ\tSN:([^\t]+).*/$1/g;
	}
	close($filedict);
}
elsif($listchrs ne "")
{
	open(my $filechrs,$listchrs);
	@chrs=<$filechrs>;
	close($filechrs);
}
else
{
	die "Invalid --ref, --dict and --listchrs options.\n$usage";
}
chomp(@chrs);
@chrs{@chrs}=(0) x @chrs;

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
my %ichrs;

# Mask variables. Arrays with booleans in the index of files that have the current chr or pos. Array-based masks
my @chrmask=(1) x $nfiles; #Initially we assume all samples are valid for chromosome and position
my @posmask=@chrmask;
my @posvalidfilesmask=@chrmask;
my @validfilesmask=@chrmask;

# Content variables
my $line;
my @results;
my @temp;

#Results init
for ($i=0; $i <= $nfiles; ++$i)
{
	for (my $j=0; $j< $nfilters; ++$j)
	{
		$results[$i][$j]=0;
	}
}

#Obtaining the first line
#$pos=getnext(@posvalidfilesmask); #Equivalent, but a little slower
for (my $nfile=0; $nfile<$nfiles; ++$nfile)
{
	if ($posvalidfilesmask[$nfile]==1)
	{
		$line=$files[$nfile]->getline();
		chomp($line);
		($cchr[$nfile],$cpos[$nfile],$cdepth[$nfile])=split("\t",$line);
	}
}

# Control variables
my $stop=0;
$i=0;
my $j=0;
my $n=0;
my $nchrs=0;
my $nchr=0;
my $chr=$chrs[$nchr];
my $pos=0;
my $ifilt=0;
my $isample=0;
my $nvalidfiles=0;
my $npassthisfilter=0;

##Main loop that gets new positions until the end of all files or an iterator limit (to avoid infinite loops with a while)
while (!($stop == $nfiles || $i>=$MAX_IT))
{
	#Check if the new positions are in this chr, if they are not, move to the next until they are
	do
	{
		#Make chr mask and check how many samples still in this chr
		#@chrmask=makemask(\@cchr,$chr);#Equivalent, but a little slower
		@chrmask=(0) x $nfiles;	
		for ($j=0; $j<$nfiles; ++$j)
		{
			if(defined $cchr[$j] and $cchr[$j] eq $chr)
			{
				$chrmask[$j]=1;
			}
		}
		
		$nchrs=sum(@chrmask);
		
		#If none, next chr
		if($nchrs==0){++$nchr};

		$nchr >= scalar @chrs and die "Unknown chromosomes present in the input files. Make sure that the input reference is the same you used to generate the bams you run with samtools depth\n";
		$chr=$chrs[$nchr];
	}
	while($nchrs ==0);

	#$pos=min(@cpos[makeslicefrommask(\@posvalidfilesmask)]); #Equivalent, but a little slower
	@temp=();
	for($j=0;$j<scalar $nfiles;++$j)
	{
		if($chrmask[$j]==1)
		{
			push(@temp,$j);
		}
	}
	$pos=min(@cpos[@temp]);
	
	#@posmask=makemask(\@cpos,$pos);#Equivalent, but a little slower
	@posmask=(0) x $nfiles;	
	for ($j=0; $j<$nfiles; ++$j)
	{
		if(defined $cpos[$j] and $cpos[$j] ==$pos)
		{
			$posmask[$j]=1;
		}
	}

	#@posvalidfilesmask=arrayfunc(\@posmask, \&vand, \@chrmask);
	#Equivalent, but this way is a little faster (although less beautiful)
	for ($j=0; $j<$nfiles; ++$j)
	{
		$posvalidfilesmask[$j]=$posmask[$j] & $chrmask[$j];
	}
	
	@validfilesmask=@posvalidfilesmask;
	$nvalidfiles=sum(@validfilesmask);

	for ($ifilt=0; $ifilt< scalar @filters; ++$ifilt) #desired depth
	{
		$npassthisfilter=0; #number of samples that pass this filter

		for($isample=0; $isample<$nfiles; ++$isample)
		{
			if($validfilesmask[$isample]==1) #Sample that is not masked out
			{
				if($cdepth[$isample]>=$filters[$ifilt]) #Depth comparison
				{	
					++$npassthisfilter;
				}
				else
				{
					$validfilesmask[$isample]=0; #this sample does not have enough coverage at this level and therefore it does not need to be checked again for this position
					--$nvalidfiles;
				}
			}
		}
		for ($j=1; $j<=$npassthisfilter; ++$j)
		{
			++$results[$j][$ifilt]; #Results for this position
		}
	}	
	
	#Loading data for the next iteration
	#$pos=getnext(@posvalidfilesmask); #Equivalent, but a little slower
	for (my $nfile=0; $nfile<$nfiles; ++$nfile)
	{
		if ($posvalidfilesmask[$nfile]==1)
		{
			$line=$files[$nfile]->getline();
			if (defined $line)
			{
				chomp($line);
				($cchr[$nfile],$cpos[$nfile],$cdepth[$nfile])=split("\t",$line);
			}
			else
			{
				($cchr[$nfile],$cpos[$nfile],$cdepth[$nfile])=(undef) x 4;
				++$stop;
			}
		}
	}
	if ($i%1000000 == 0)
	{
		print("Iteraction number $i, CHR $chr, POS $pos\n");
	}
	++$i;
}

# Print results
# #############

my @header;

##MATRIX output
###############
push(@header,">=$_") foreach @filters;

print(join(",","Sharedby",@header),"\n");
for ($i=1; $i <= $nfiles; ++$i)
{
	print(join(",",$i,@{$results[$i]}),"\n");
}

##Tidy output (i.e., R-style)
#############################

open(my $ofilehand,">$output_prefix.csv");
@header=qw(bp mindepth nsamples);

print($ofilehand join(",",@header),"\n");
for ($i=1; $i<= $nfiles; ++$i)
{
	for($j=0; $j< $nfilters; ++$j){
		print($ofilehand join(",",$results[$i][$j],$filters[$j],$i),"\n");
	}
}
close($ofilehand);

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

sub makeslicefrommask
{
	my @out;

	for(my $i=0;$i<scalar @{$_[0]};++$i)
	{
		if($_[0]->[$i]==1)
		{
			push(@out,$i);
		}
	}

	return @out;
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

sub getnext
{
	my @filemask=@_;

	for (my $nfile=0; $nfile<$nfiles; ++$nfile)
	{
		if ($filemask[$nfile]==1)
		{
			$line=$files[$nfile]->getline();
			chomp($line);
			($cchr[$nfile],$cpos[$nfile],$cdepth[$nfile])=split("\t",$line);
			#print("DEBUG: chr=$cchr[$nfile], pos=$cpos[$nfile], depth=$cdepth[$nfile]\n");
		}
	}
	return min(@cpos);
}

##By amon at https://stackoverflow.com/questions/18103501/prompting-multiple-questions-to-user-yes-no-file-name-input
sub prompt {
	my ($query) = @_; # take a prompt string as argument
 	local $| = 1; # activate autoflush to immediately show the prompt
	print $query;
	chomp(my $answer = <STDIN>);
	return $answer;
}

##By Xtof at https://stackoverflow.com/questions/18103501/prompting-multiple-questions-to-user-yes-no-file-name-input
sub prompt_yn {
	my ($query, $default) = @_; 
	my $default_yes = lc $default eq 'y';
	my $yn = $default_yes ? "[Y/n]" : "[y/N]";
	my $answer = lc prompt("$query $yn: "); 
	return $default_yes ? ! ($answer =~ /^n/) : $answer =~ /^y/;
}
