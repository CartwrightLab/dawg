#! /usr/bin/perl -w
# Copyright (2005) Reed A. Cartwright.  All rights reserved.
#
# converts nexus sequences to fasta
#
# usage: perl nexus2fasta.pl < infile > outfile
#
# Distributed under the same license as DAWG
#


my $state = 0;

my %seqs = ();

local $/;
my $text = <>;

my ($data) = $text =~ /begin\s+data;.+?matrix\s*(.*?);\s*end;/is;
my @lines = split(/\n/, $data);

foreach(@lines)
{
	s/^\s+//;
	s/\s+$//;
	next unless(/\w/);
	my @sec = split(/\s+/, $_);
	my $name = shift(@sec);
	$seqs{$name} |= '';
	$seqs{$name} .= join('', @sec);
}

print ">$_\n$seqs{$_}\n\n" foreach(sort(keys(%seqs)));

