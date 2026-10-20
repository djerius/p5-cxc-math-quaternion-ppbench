use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);
my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);
sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 1; (time-$t0)/$n }
for my $N (1e6, 3e6, 1e7) { my $v = random(3,$N)-0.5; printf "%9.0e new out %6.1f ns\n", $N, 1e9*bench(sub { rotq($q,$v) })/$N }
