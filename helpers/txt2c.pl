#! /usr/bin/perl -w

use strict;

while(<>)
{
	chomp;
	s/\s+$//;
	s/\\/\\\\/;
	s/"/\\"/;
	print qq("$_\\n" \\\n);
}