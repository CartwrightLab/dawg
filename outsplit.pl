#! /usr/bin/perl -w
# Copyright (2005) Reed A. Cartwright.  All rights reserved.
#
# outsplit.pl is used to extract sequences from Fasta and Phylip files
# 
# usage: perl -w outsplit.pl <file> <id>
#
# if <id> is "all" it it creates a directory and
# creates a new file for each alignment
#
# Distributed under the same license as DAWG
#

use strict;
use File::Basename;
use File::Path;
use File::Spec::Functions;

my ($file, $id) = @ARGV;

local $/;

open( FILE, $file) or die("Error opening $file.");

my $text = <FILE>;

close(FILE);

my @blocks = split(/\[DataSet \d+\].*\n/, $text);

if($id ne 'all')
{
	print $blocks[$id];	
}
else
{
	my ($name,$dir,$ext) = fileparse($file, qr{\..*});
	my $outdir = catdir($dir, $name);
	print "Creating directory $outdir\n";
	mkpath($outdir) unless (-d $outdir);
	chdir($outdir) or die("Unable to change directory.");
	
	foreach my $i (1..$#blocks)
	{
		my $out = "${name}_$i$ext";
		print "Creating file $out\n";
		open(OUT, ">$out") or die("Unable to open file.");
		print OUT $blocks[0];
		print OUT $blocks[$i];
		close(OUT);
	}
}
