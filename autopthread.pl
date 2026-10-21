use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# rotq (preallocated output) at sizes below PDL's default autopthread size
# threshold, with the default threshold and with the threshold set to 0, so
# that PDL splits the work across threads at every size.  Both are timed back
# to back in the same process.

my $q = pdl(0.1, -0.3, 0.5, 0.8); $q /= sqrt(($q**2)->sum);
sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

my $default = PDL::get_autopthread_size();
printf "PDL %s, autopthread targ = %d, default size threshold = %s M elements\n\n",
    $PDL::VERSION, PDL::get_autopthread_targ(), $default;
printf "%9s %14s %14s\n", 'N', 'default size', 'size 0';
for my $N (1e3, 1e4, 1e5, 3e5) {
    my $v = random(3, $N) - 0.5; my $r = zeroes(3, $N);
    PDL::set_autopthread_size($default);
    my $t1 = 1e9 * bench(sub { rotq($q, $v, $r) }) / $N;
    PDL::set_autopthread_size(0);
    my $t2 = 1e9 * bench(sub { rotq($q, $v, $r) }) / $N;
    printf "%9.0e %14.1f %14.1f\n", $N, $t1, $t2;
}
