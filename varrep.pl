#! /usr/bin/perl -w -ibak
# Copyright (2005) Reed A. Cartwright.  All rights reserved.
#
# varrep.pl is used to substitute varables in files
#
# usage: perl varrep.pl file
#
# right now the only variable supported is #NUM#
# useful for modifying the code blocks of Nexus files produced by Dawg
#
# Distributed under the same license as DAWG
#


use strict;

my %vars = (NUM => 0);

while(<>)
{
	$vars{NUM} = $1 if(/\[DataSet (\d+)\]/);
	s/#(\w+)#/$vars{$1}/ge;
	print $_;
}