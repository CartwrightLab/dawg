use strict;
use warnings;

my @n = qw(T C A G);
my %n = (T => 0, C => 1, A => 2, G => 3);
my @e = ([-1,0,1,2], [-1,-1,3,4], [-1,-1,-1,5], [-1,-1,-1,-1]);

my @d;

my @row = ();
foreach my $i(0..62) {
	my @a = @n[int($i/16)%4, int($i/4) % 4, $i % 4];
	foreach my $j(($i+1)..63) {
		my @b = @n[int($j/16)%4, int($j/4) % 4, $j % 4];
		my $x = 0;
		++$x if($a[0] ne $b[0]);
		++$x if($a[1] ne $b[1]);
		++$x if($a[2] ne $b[2]);
		if($x > 1) {
			push(@row, -1);
		} elsif($a[0] ne $b[0]) {
			push(@row, 0+$e[$n{$a[0]}][$n{$b[0]}]);
		} elsif($a[1] ne $b[1]) {
			push(@row, 8+$e[$n{$a[1]}][$n{$b[1]}]);
		} else {
			push(@row, 16+$e[$n{$a[2]}][$n{$b[2]}]);
		}
	}
}

@row = map { sprintf("% 2s", $_) } @row;

my @rrow = ();
my $y = int($#row / 24);
foreach my $r(0..$y) {
	my $a = ($r*24);
	my $b = (($r+1)*24-1);
	$b = @row-1 if($b >= @row);
	push(@rrow, join(",", @row[$a..$b]));
}

print join(",\n", @rrow) . "\n";

