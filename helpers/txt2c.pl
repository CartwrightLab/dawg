#! /usr/bin/perl -w

use strict;

while(<>)
{
	chomp;
	s/\s+$//;
	s/\\/\\\\/g;
	s/"/\\"/g;
	print qq("$_\\n" \\\n);
}
print qq(\n);
