use strict;
use warnings;

my @n = qw(T C A G);

my @d;

foreach my $i(0..63) {
	my @a = @n[int($i/16), int($i/4) % 4, $i % 4];
	my @row = ();
	foreach my $j(0..63) {
		my @b = @n[int($j/16), int($j/4) % 4, $j % 4];
		my $x = 0;
		++$x if($a[0] ne $b[0]);
		++$x if($a[1] ne $b[1]);
		++$x if($a[2] ne $b[2]);
		push(@row, $x);
	}
	push(@d, join(",", @row));
}

print join(",\n", @d);
print "\n";

