use strict;
use warnings;

my @a = qw(A C D E F G H I K L M N P Q R S T V W Y);

my @b = map {(ord($_) & 95 )-ord('@')} @a;

my @q = (20) x 32;

my $i = 0;
$q[$_] = $i++ foreach @b;

@q = map { sprintf("%2d", $_) } @q;

print join(",", @q[0..15]) . ",\n" . join(",", @q[16..31]) . "\n";

