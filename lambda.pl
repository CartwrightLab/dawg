#! /usr/bin/perl -w

# lambda.pl v1.0 - An estimator for DAWG's lambda
# Copyright (2004) Reed A. Cartwright.  All rights reserved.
#
# Usage: perl lambda.pl [treefile] [fastafile]
#
# Distributed under the same license as DAWG
#

use strict;

my $gettree = 1;
my $tree = '';
my %seqs = ();
my $seqname = '';

# process STDIN line by line
while(<>)
{
	# remove newline and flanking whitespace
	chomp;   
	s/^\s+//;
	s/\s+$//;
	#check to see if the programs should shift to reading sequences		
	$gettree = 0 if($gettree && m/^>/);
	if($gettree)
	{
		$tree .= $_;
	}
	#check to see if the line is a fasta sequence name
	elsif(m/^>\s*(.+)\s*/)
	{
		$seqname = $1;
	}
	#append the line to the sequence
	else
	{
		$seqs{$seqname} .= $_;
	}
}

#calculate the total length of the tree by extracting all the branch lengths
my $totlen = 0;
$totlen += $1 while($tree =~ m/:([\d\.]+)/g);

#convert sequences to a table
my @table = values(%seqs);

my %gaps = ();

my $avgL = 0;
#find unique gaps and count the ungapped length of the sequences
foreach(@table)
{
	$gaps{$-[0], $+[0]} = 1 while($_ =~ m/-+/g);
	$avgL += tr/ACGTacgt//;
}
#find average length
$avgL /= @table;

#count gaps
my $numgaps = keys(%gaps);

#calculate lambda estimate

my $lambda = $numgaps/$avgL/$totlen;
print "Total Len is $totlen.\n";
print "Number of gaps is $numgaps.\n";
print "Average Length is $avgL.\n";
print "Lambda is $lambda.\n";