use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use PDL::LinearAlgebra::Real qw(mmult gemm); use QBench; use Time::HiRes qw(time);

# PDL::LinearAlgebra (BLAS gemm via mmult / gemm) vs. the PP quaternion kernel.
# Run with and without OPENBLAS_NUM_THREADS=1.
#
# Loading PDL::LinearAlgebra::Real replaces the overload of PDL's `x` operator
# with its BLAS-based mmult, so `x` below is BLAS; PDL's own matrix multiply
# is called explicitly as PDL::matmult.

my $q = pdl(0.1, -0.3, 0.5, 0.8); $q /= sqrt(($q**2)->sum);
my $R  = matrices($q);                     # (3,3): $R x ... convention as in bench.pl
my $Rt = $R->transpose->copy;

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

# rotation matrices from quaternions (4,N...) -> (3,3,N...), with ordinary PDL operations
sub matrices {
    my ($qs) = @_;
    my ($x, $y, $z, $w) = map { $qs->slice("($_)") } 0..3;
    my @d = $qs->dims; shift @d;
    my $m = zeroes($qs->type, 3, 3, @d);
    $m->slice('(0),(0)') .= 1-2*($y*$y+$z*$z); $m->slice('(1),(0)') .= 2*($x*$y-$z*$w);   $m->slice('(2),(0)') .= 2*($x*$z+$y*$w);
    $m->slice('(0),(1)') .= 2*($x*$y+$z*$w);   $m->slice('(1),(1)') .= 1-2*($x*$x+$z*$z); $m->slice('(2),(1)') .= 2*($y*$z-$x*$w);
    $m->slice('(0),(2)') .= 2*($x*$z-$y*$w);   $m->slice('(1),(2)') .= 2*($y*$z+$x*$w);   $m->slice('(2),(2)') .= 1-2*($x*$x+$y*$y);
    return $m;
}

# correctness check against rotq
{
    my $v = random(3, 1000) - 0.5;
    my $ref = rotq($q, $v);
    my $c = zeroes(3, 1000); gemm($R, 1, 0, $v, 1, 0, $c);
    my $qs = random(4, 1000) - 0.5; $qs /= sqrt(($qs**2)->sumover)->dummy(0);
    printf "max |diff| vs rotq: matmult %.1e  x %.1e  mmult %.1e  gemm %.1e  broadcast mmult %.1e\n",
        (abs(PDL::matmult($v, $Rt) - $ref))->max, (abs($v x $Rt - $ref))->max, (abs(mmult($v, $Rt) - $ref))->max, (abs($c - $ref))->max,
        (abs(mmult($v->dummy(1), matrices($qs)->xchg(0,1))->slice(':,(0)') - rotq($qs, $v)))->max;
}
printf "PDL %s, PDL::LinearAlgebra %s, autopthread targ = %d, OPENBLAS_NUM_THREADS=%s\n\n",
    $PDL::VERSION, $PDL::LinearAlgebra::Real::VERSION // $PDL::LinearAlgebra::VERSION // '?',
    PDL::get_autopthread_targ(), $ENV{OPENBLAS_NUM_THREADS} // 'default';

print "(1) one quaternion, N vectors: ns per vector\n";
printf "%9s %9s %11s %9s %11s %9s %12s\n", 'N', 'rotq', 'rotq pre', 'matmult', 'x (=mmult)', 'mmult', 'gemm pre';
for my $N (1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7) {
    my $v = random(3, $N) - 0.5;
    my ($r, $c) = (zeroes(3, $N), zeroes(3, $N));
    printf "%9.0e %9.1f %11.1f %9.1f %11.1f %9.1f %12.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($q, $v) },
        sub { rotq($q, $v, $r) },
        sub { PDL::matmult($v, $Rt) },
        sub { $v x $Rt },
        sub { mmult($v, $Rt) },
        sub { gemm($R, 1, 0, $v, 1, 0, $c) };
}

print "\n(2) N quaternions, one vector each: ns per vector\n";
printf "%9s %9s | %16s %16s\n", 'N', 'rotq', 'mmult (mats pre)', 'build + mmult';
for my $N (1e2, 1e4, 1e6) {
    my $qs = random(4, $N) - 0.5; $qs /= sqrt(($qs**2)->sumover)->dummy(0);
    my $v = random(3, $N) - 0.5;
    my $mt = matrices($qs)->xchg(0,1)->copy;          # (3,3,N), transposed, built once
    printf "%9.0e %9.1f | %16.1f %16.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($qs, $v) },
        sub { mmult($v->dummy(1), $mt) },
        sub { mmult($v->dummy(1), matrices($qs)->xchg(0,1)) };
}
