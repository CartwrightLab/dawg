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

#calculate distribution of gaps sizes
my %gapsizes = ();
my $maxgap = 0;
my $suml = 0;
foreach(keys(%gaps))
{
	my ($f, $l) = split(/$;/);
	$l -= $f;
	$gapsizes{$l} ||=0;
	$gapsizes{$l}++;
	$maxgap = $l if($l > $maxgap);
	$suml += $l;
}

#use liklihood to estimate p and q
my $avgG = $suml/$numgaps-1;
my %qhat = ();
my $rlog = 0;
foreach my $r(1..$maxgap)
{	
	#estimate q | r
	my $q = $avgG/($avgG+$r);
	#calculate LL(q,r)
	my $LL = $numgaps*($r*log(1-$q)-$rlog-log($q));
	while(my($g,$n) = each(%gapsizes))
	{
		$LL += $n*$g*log($q);
		$LL += $n*log($_) foreach($g..$g+$r-2);
	}
	#calculate BIC to test whether r=1 is a better explaination
	my $bic = -2*$LL+log($numgaps)*($r==1 ? 1 : 2);
	#store data
	$qhat{$r} = [$q, $LL, $bic];
	$rlog += log($r);
}
#find maximimum liklihood and minimium BIC
my $r_bic = 1;
my $r_ll = 1;
foreach my $r(1..$maxgap)
{
	$r_bic = $r if($qhat{$r}->[2] < $qhat{$r_bic}->[2]);
	$r_ll = $r if($qhat{$r}->[1] > $qhat{$r_ll}->[1]);
}

#output
print "Total Len is $totlen.\n";
print "Number of gaps is $numgaps.\n";
print "Average Length is $avgL.\n";
print "\nLambda Estimate is $lambda.\n";
print "\nGap Size Distribution:\n";
print join("\t", $_, $gapsizes{$_} || 0, ($gapsizes{$_} || 0)/$numgaps), "\n" foreach(1..$maxgap);
print "\nEstimates of r and q based on Maximum Liklihood\n";
print "\tr Estimate is $r_ll.  (BIC $qhat{$r_ll}->[2])\n";
print "\tq Estimate is ${qhat{$r_ll}->[0]}.  (LH ${qhat{$r_ll}->[2]})\n";
print "\nEstimates of r and q based on Minimum BIC\n";
print "\tr Estimate is $r_bic.  (BIC $qhat{$r_bic}->[2])\n";
print "\tq Estimate is $qhat{$r_bic}->[0].  (LH $qhat{$r_bic}->[2])\n";
