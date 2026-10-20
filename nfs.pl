use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# Cost of allocating an output with zeroes vs. new_from_specification (no fill),
# and of passing either, freshly allocated on every call, to rotq / rotq3.

my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);
sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }
sub minflt { open my $f, '<', '/proc/self/stat' or die; (split ' ', scalar <$f>)[9] }
sub nfs { my $p = PDL->new_from_specification(double, @_); $p->make_physical; $p }

printf "PDL %s, MALLOC_MMAP_THRESHOLD_=%s\n", $PDL::VERSION, $ENV{MALLOC_MMAP_THRESHOLD_} // 'default';
{
    my $p = nfs(3, 1e6);
    printf "new_from_specification: allocated=%d  (contents not guaranteed; max |x| here = %g)\n\n",
        $p->allocated, $p->abs->max;
}

print "allocation alone: ns per vector (3 doubles), and minor page faults per call\n";
printf "%9s %9s %9s %9s %9s\n", 'N', 'zeroes', 'nfs', 'flt zero', 'flt nfs';
for my $N (1e3, 1e5, 1e6, 3e6, 1e7) {
    my @t = map { 1e9 * bench($_) / $N } sub { zeroes(3,$N) }, sub { nfs(3,$N) };
    my @f = map { my $c = $_; $c->(); my $f0 = minflt(); $c->() for 1..10; (minflt()-$f0)/10 }
        sub { zeroes(3,$N) }, sub { nfs(3,$N) };
    printf "%9.0e %9.2f %9.2f %9.0f %9.0f\n", $N, @t, @f;
}

print "\nrotating with a freshly allocated output on every call: ns per vector\n";
printf "%9s | %9s %11s %9s %9s | %9s %11s %9s %9s\n", 'N',
    'rotq', 'rotq zeroes', 'rotq nfs', 'rotq pre', 'rotq3', 'rotq3 zeroes', 'rotq3 nfs', 'rotq3 pre';
for my $N (1e3, 1e5, 1e6, 3e6, 1e7) {
    my $v = random(3,$N) - 0.5; my @c = map { random($N) - 0.5 } 1..3;
    my $r = zeroes(3,$N); my @o = map { zeroes($N) } 1..3;
    printf "%9.0e | %9.1f %11.1f %9.1f %9.1f | %9.1f %11.1f %9.1f %9.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($q, $v) },
        sub { rotq($q, $v, zeroes(3,$N)) },
        sub { rotq($q, $v, nfs(3,$N)) },
        sub { rotq($q, $v, $r) },
        sub { my @r = rotq3($q, @c) },
        sub { rotq3($q, @c, map { zeroes($N) } 1..3) },
        sub { rotq3($q, @c, map { nfs($N) } 1..3) },
        sub { rotq3($q, @c, @o) };
}
