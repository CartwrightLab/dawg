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
my $avgG = $suml/$numgaps-1.0;

my %model = ();

#geometric model
#MLE of q
my $q = $avgG/($avgG+1.0);
#calculate LL(q)
my $LL = $numgaps*(log(1.0-$q)-log($q))+$suml * log($q);
#calculate BIC 
my $bic = -2*$LL+log($numgaps)*1.0;
#calculate xsq
my $xsq = 0.0;
my $df = $maxgap-2;
while(my($g,$n) = each(%gapsizes))
{
	my $e = $numgaps*((1.0-$q) * $q**($g-1.0));
	$xsq += ($n-$e)**2/$e;
}


$model{Geometric} = {BIC => $bic, LL => $LL, DF => $df, XSQ => $xsq, Params => {q => $q}};

my %rmodel = ();
my $rlog = 0;
foreach my $r(2..$maxgap)
{	
	#estimate q | r
	my $q = $avgG/($avgG+$r);
	#calculate LL(q,r)
	my $LL = $numgaps*($r*log(1-$q)-$rlog-log($q))+$suml * log($q);
	my $xsq = 0.0;
	my $df = $maxgap-3;
	while(my($g,$n) = each(%gapsizes))
	{
		my $sl = 0.0;
		$sl += log($_) foreach($g..$g+$r-2);
		my $e = $numgaps*( (1.0-$q)**$r * $q**($g-1) /exp($rlog)*exp($sl) );
		$xsq += ($n-$e)**2/$e;
		$LL += $n*$sl;
	}
	#calculate BIC
	my $bic = -2*$LL+log($numgaps)*2.0;
	#store data
	$rmodel{$r} = {BIC => $bic, LL => $LL,  DF => $df, XSQ => $xsq, Params => {r => $r, q => $q }};
	$rlog += log($r);
}
#find maximimum liklihood
my $r_ll = 2;
foreach my $r(3..$maxgap)
{
	$r_ll = $r if($rmodel{$r}{LL} > $rmodel{$r_ll}{LL});
}
$model{NegBin} = $rmodel{$r_ll};

# truncated power-law
my @gs = (sort(keys(%gapsizes)))[0..4];

my $sx = 0.0;
my $sx2 = 0.0;
my $sxy = 0.0;
my $sy = 0.0;
my $n = @gs;
foreach(@gs)
{
	my $x = log($_);
	my $y = log($gapsizes{$_}/$numgaps);
	$sx += $x;
	$sy += $y;
	$sxy += $x*$y;
	$sx2 += $x*$x;
}
my $a = ($n*$sxy-$sx*$sy)/($n*$sx2-$sx*$sx);
my $b = ($sy-$a*$sx)/$n;
$model{PowerLaw} = {BIC => $bic, LL => $LL, DF => $df, XSQ => $xsq, Params => {a => -$a, b => exp($b)}};

#output
print "Total Len is $totlen.\n";
print "Number of gaps is $numgaps.\n";
print "Average Length is $avgL.\n";
print "\nLambda Estimate is $lambda.\n";
print "\nAverage Gap Size is ", $avgG + 1.0, ".\n";
print "\nGap Size Distribution:\n";
print join("\t", $_, $gapsizes{$_} || 0, ($gapsizes{$_} || 0)/$numgaps), "\n" foreach(1..$maxgap);
print "\nEstimated Models\n";
foreach my $k (sort(keys(%model)))
{
	print "$k:\n\t";
	my @par = ();
	while(my($p,$v) = each(%{$model{$k}{Params}}))
	{
		push(@par, "$p = $v");
	}
	print join(', ', @par), "\n";
	print "\tLogLik = $model{$k}{LL}\n";
	print "\tBIC    = $model{$k}{BIC}\n";
	print "\tXSQ    = $model{$k}{XSQ} ($model{$k}{DF} df)\n";
}

