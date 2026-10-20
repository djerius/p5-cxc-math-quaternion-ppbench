use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);
my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);
sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 1; (time-$t0)/$n }
printf "autopthread size threshold: %s M elements\n", PDL::get_autopthread_size();
printf "%9s %10s %10s %12s\n", 'N', 'new out', 'prealloc', '1 thread pre';
for my $N (1e5, 1e6, 3e6, 1e7, 3e7) {
    my $v = random(3,$N) - 0.5; my $out = zeroes(3,$N);
    my @r = map { 1e9*bench($_)/$N } sub { rotq($q,$v) }, sub { rotq($q,$v,$out) };
    PDL::set_autopthread_targ(1); push @r, 1e9*bench(sub { rotq($q,$v,$out) })/$N; PDL::set_autopthread_targ(8);
    printf "%9.0e %10.1f %10.1f %12.1f\n", $N, @r;
}
