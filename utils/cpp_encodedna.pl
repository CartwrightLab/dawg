use strict;
use warnings;

my @a = qw(A C G T);
my @b = map {ord($_)-ord('0')} @a;
my @q = (63) x 80;

my $j = 0;

my @x = (63) x 80;
foreach(@a) {
	$x[ord(lc($_))-ord('0')] = $j;
	$x[ord($_)-ord('0')] = $j;
	$j++;
}
$x[ord('u')-ord('0')] = 3;
$x[ord('U')-ord('0')] = 3;

@x = map { sprintf("% 2s", $_) } @x;
print join(",", @x[ 0..19]) . ",\n" . 
	join(",", @x[20..39]) . ",\n" .
	join(",", @x[40..59]) . ",\n" .
	join(",", @x[60..79]) . "\n"
  ;


