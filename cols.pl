use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# Vectors supplied as three independent columns $x, $y, $z, each of shape (N).
# Compare (1) building a (3,N) ndarray from them, and (2) rotating them,
# either via a (3,N) ndarray and rotq, or directly with the column kernel rotq3.

my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

# ways to build a physical (3,N) ndarray from the columns
my %combine = (
    'cat+xchg+copy' => sub { cat(@_)->xchg(0,1)->copy },
    'glue'          => sub { $_[0]->dummy(0)->glue(0, $_[1]->dummy(0), $_[2]->dummy(0)) },
);
# fill a preallocated (3,N) ndarray
sub fill { my ($v, @c) = @_; $v->slice("($_)") .= $c[$_] for 0..2; $v }

# correctness check
{
    my ($x,$y,$z) = map { random(1000) - 0.5 } 1..3;
    my $ref = rotq($q, cat($x,$y,$z)->xchg(0,1)->copy);
    my @r3 = rotq3($q, $x, $y, $z);
    my $v = $combine{glue}->($x,$y,$z);
    printf "max |diff|: rotq3 %.1e  glue %.1e  xchg-view %.1e  fill %.1e\n",
        (abs(cat(@r3)->xchg(0,1) - $ref))->max,
        (abs(rotq($q,$v) - $ref))->max,
        (abs(rotq($q, cat($x,$y,$z)->xchg(0,1)) - $ref))->max,
        (abs(rotq($q, fill(zeroes(3,1000),$x,$y,$z)) - $ref))->max;
}
printf "PDL %s, autopthread targ = %d, MALLOC_MMAP_THRESHOLD_=%s\n\n",
    $PDL::VERSION, PDL::get_autopthread_targ(), $ENV{MALLOC_MMAP_THRESHOLD_} // 'default';

my @Ns = (1e2, 1e3, 1e4, 1e5, 1e6, 1e7);

print "(1) building a (3,N) ndarray from three (N) columns: ns per vector\n";
printf "%9s %14s %9s %13s\n", 'N', 'cat+xchg+copy', 'glue', 'fill prealloc';
for my $N (@Ns) {
    my @c = map { random($N) - 0.5 } 1..3; my $v = zeroes(3,$N);
    printf "%9.0e %14.1f %9.1f %13.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { $combine{'cat+xchg+copy'}->(@c) }, sub { $combine{glue}->(@c) },
        sub { fill($v, @c) };
}

print "\n(2) rotating: ns per vector\n";
printf "%9s %9s %11s | %11s %14s %9s %13s | %9s %13s\n", 'N',
    'rotq(3,N)', 'rotq(N,3)^T', 'cat+xchg', 'cat+xchg+copy', 'glue', 'fill prealloc', 'rotq3', 'rotq3 prealloc';
printf "%9s %21s | %51s | %23s\n", '', '--- no build ---', '------- build (3,N), then rotq -------', '---- columns direct ----';
for my $N (@Ns) {
    my @c = map { random($N) - 0.5 } 1..3;
    my $v0 = cat(@c)->xchg(0,1)->copy;            # already in (3,N) layout
    my $vt = cat(@c);                             # (N,3): columns stacked, a 2-D column store
    my ($vp, $rp) = (zeroes(3,$N), zeroes(3,$N));  # preallocated (3,N) in/out
    my @o = map { zeroes($N) } 1..3;              # preallocated output columns
    printf "%9.0e %9.1f %11.1f | %11.1f %14.1f %9.1f %13.1f | %9.1f %13.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($q, $v0) },
        sub { rotq($q, $vt->xchg(0,1)) },        # transposed view of (N,3), no copy
        sub { rotq($q, cat(@c)->xchg(0,1)) },      # transposed view, no copy
        sub { rotq($q, $combine{'cat+xchg+copy'}->(@c)) },
        sub { rotq($q, $combine{glue}->(@c)) },
        sub { rotq($q, fill($vp, @c), $rp) },
        sub { my @r = rotq3($q, @c) },
        sub { rotq3($q, @c, @o) };
}
