use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL;
use PDL::Transform;
use QBench;
use Time::HiRes qw(time);

# unit quaternion, scalar-last (x,y,z,w)
my $q = pdl( 0.1, -0.3, 0.5, 0.8 );
$q /= sqrt( ( $q**2 )->sum );
my ( $x, $y, $z, $w ) = $q->list;
my $R = pdl(
    [
        1 - 2 * ( $y * $y + $z * $z ),
        2 * ( $x * $y - $z * $w ),
        2 * ( $x * $z + $y * $w )
    ],
    [
        2 * ( $x * $y + $z * $w ),
        1 - 2 * ( $x * $x + $z * $z ),
        2 * ( $y * $z - $x * $w )
    ],
    [
        2 * ( $x * $z - $y * $w ),
        2 * ( $y * $z + $x * $w ),
        1 - 2 * ( $x * $x + $y * $y )
    ]
);
my $t = t_linear( matrix => $R->transpose );    # t_linear computes v x M

# plain PDL ops version of rotq (temporaries at every step), q may be (4,N)
sub rotq_ops {
    my ( $q, $v ) = @_;
    my $u = $q->slice('0:2');
    my $w = $q->slice('(3)');
    my $t = 2 * $u->crossp($v);
    return $v + $w->dummy(0) * $t + $u->crossp($t);
}

sub bench {
    my ($code) = @_;
    $code->();    # warm-up
    my ( $n, $t0 ) = ( 0, time );
    do { $code->(); $n++ } while time - $t0 < 0.5;
    return ( time - $t0 ) / $n;
}

# correctness check
{
    my $v   = random( 3, 1000 ) - 0.5;
    my $ref = rotq( $q, $v );
    printf
      "max |diff| vs rotq: rotm %.1e  matmult %.1e  t_linear %.1e  ops %.1e\n",
      ( abs( rotm( $R, $v ) - $ref ) )->max,
      ( abs( $v x $R->transpose - $ref ) )->max,
      ( abs( $t->apply($v) - $ref ) )->max,
      ( abs( rotq_ops( $q, $v ) - $ref ) )->max;
}

printf "PDL %s, autopthread targ = %d\n\n", $PDL::VERSION,
  PDL::get_autopthread_targ();
print "one quaternion, N vectors: ns per vector\n";
printf "%9s %9s %9s %9s %9s %9s\n", 'N', 'rotq(PP)', 'rotm(PP)', 'matmult',
  't_linear', 'PDL ops';
for my $N ( 1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7 ) {
    my $v = random( 3, $N ) - 0.5;
    my @r = map { 1e9 * bench($_) / $N } sub { rotq( $q, $v ) },
      sub { rotm( $R, $v ) }, sub { $v x $R->transpose },
      sub { $t->apply($v) }, sub { rotq_ops( $q, $v ) };
    printf "%9.0e" . " %9.1f" x 5 . "\n", $N, @r;
}

print "\nN quaternions, N vectors (one each): ns per vector\n";
printf "%9s %9s %9s\n", 'N', 'rotq(PP)', 'PDL ops';
for my $N ( 1e3, 1e5, 1e7 ) {
    my $qs = random( 4, $N ) - 0.5;
    $qs /= sqrt( ( $qs**2 )->sumover )->dummy(0);
    my $v = random( 3, $N ) - 0.5;
    printf "%9.0e %9.1f %9.1f\n", $N,
      map { 1e9 * bench($_) / $N } sub { rotq( $qs, $v ) },
      sub { rotq_ops( $qs, $v ) };
}
