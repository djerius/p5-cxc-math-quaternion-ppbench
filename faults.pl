use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench;
# minor page faults per rotq call (field 10 of /proc/self/stat)
sub minflt { open my $f, '<', '/proc/self/stat' or die; my @s = split ' ', scalar <$f>; $s[9] }
my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);
printf "%9s %14s %14s %10s\n", 'N', 'faults/call', 'out pages', 'prealloc';
for my $N (1e6, 3e6, 1e7) {
    my $v = random(3,$N) - 0.5; my $out = zeroes(3,$N);
    rotq($q,$v); rotq($q,$v,$out);
    my $f0 = minflt(); rotq($q,$v) for 1..20; my $new = (minflt()-$f0)/20;
    $f0 = minflt(); rotq($q,$v,$out) for 1..20; my $pre = (minflt()-$f0)/20;
    printf "%9.0e %14.0f %14.0f %10.0f\n", $N, $new, 3*8*$N/4096, $pre;
}
