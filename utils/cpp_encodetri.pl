use strict;
use warnings;

my @a = split(//, 'ABCDEFGHIJ_:KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789');

my @b = map {ord($_)-ord('0');} @a;

my @q = (10) x 80;

my $i = 0;
$q[$_] = $i++ foreach @b;

@q = map { sprintf("%2d", $_) } @q;

print join(",", @q[ 0..15]) . ",\n" . 
      join(",", @q[16..31]) . ",\n" .
      join(",", @q[32..47]) . ",\n" .
      join(",", @q[48..63]) . ",\n" .
      join(",", @q[64..79]) . "\n"
      ;

