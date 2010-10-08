use strict;
use warnings;

my @a = split(//, 'ABCDEFGHIJ_:KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789');

my @b = map {ord($_)-ord('0');} @a;

my @q = (10) x 80;

my $i = 0;
$q[$_] = $i++ foreach @b;

@q = map { sprintf("%2d", $_) } @q;

print "\n// char -> triplet\n";
print join(",", @q[ 0..15]) . ",\n" . 
      join(",", @q[16..31]) . ",\n" .
      join(",", @q[32..47]) . ",\n" .
      join(",", @q[48..63]) . ",\n" .
      join(",", @q[64..79]) . "\n"
      ;

my @nord = ('T', 'C', 'A', 'G'); 	  
my @codes = (
	# 1 The Standard Code
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    # 2 The Vertebrate Mitochondrial Code
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    # 3 The Yeast Mitochondrial Code
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    # 4 The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	# 5 The Invertebrate Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
	# 6 The Ciliate, Dasycladacean and Hexamita Nuclear Code
	"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	"",
	"",
	# 9 The Echinoderm and Flatworm Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	#10 The Euplotid Nuclear Code
	"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	#11 The Bacterial, Archaeal, and Plant Plastid Code
	"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	#12 The Alternative Yeast Nuclear Code
	"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	#13 The Ascidian Mitochondrial Code
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
	#14 The Alternative Flatworm Mitochondrial Code
	"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	#15 The Blepharisma Nuclear Code
	"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	#16 The Chlorophycean Mitochondrial Code
	"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	"",
	"",
	"",
	"",
	#21 The Trematode Mitochondrial Code 
	"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
	#22 The Scenedesmus obliquus mitochondrial Code 
	"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
	#23 The Thraustochytrium Mitochondrial Code
	"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
);

@a = split(//, 'ABCDEFGHIJ_:KLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789!-');


print("\n// base -> char\n");
foreach(@codes) {
	my @c = split(//);
	my $j = 0;
	my $i = 0;
	my @x = (64) x 64;
	foreach(@c) {
		$x[$j++] = $i if($_ ne '*');
		$i++;
	}
	$x[63] = 65;
	print '"' . join('', @a[@x]) . "\"\n";	
}

print("\n// char -> base\n");
foreach(@codes) {
	my @c = split(//);
	my $j = 0;
	my $i = 0;
	my @x = (63) x 80;
	foreach(@c) {
		$x[ord($a[$i])-ord('0')] = $j++ if($_ ne '*');
		$i++;
	}
	@x = map { sprintf("% 2s", $_) } @x;
	print join(",", @x[ 0..19]) . ",\n" . 
		join(",", @x[20..39]) . ",\n" .
		join(",", @x[40..59]) . ",\n" .
		join(",", @x[60..79]) . "\n"
      ;
}

print("\n// removing stops\n");
foreach(@codes) {
	my @x = ();
	while(/([*])/g) {
		push(@x, pos()-length($1));
	}
	@x = reverse(@x);
	push(@x, 0) while(@x < 5);
	@x = map { sprintf("% 2s", $_) } @x;
	print "\t\t" . join(',', @x) . ",\n";
}



