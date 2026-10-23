use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use Scalar::Util qw( refaddr );

# Checks the generic dispatch engine (qdispatch.h; qd_rotate_vec, qd_mult,
# qd_pow, qd_to_matrix) against the PP kernels and against identities, for
# every combination of input forms and output modes.

my ( $n, $fail ) = ( 0, 0 );
sub ok { my ( $c, $name ) = @_; $n++; return 1 if $c; $fail++; print "FAIL: $name\n"; 0 }

# any result form -> ndarray (values along dim 0, like the PP kernels)
sub to_pdl {
    my ( $x, $keys ) = @_;
    return $x if ref $x && eval { $x->isa('PDL') };
    return cat( map { pdl( $x->{$_} ) } @$keys )->xchg( 0, 1 ) if ref $x eq 'HASH';
    return pdl($x);
}
sub same { my ( $a, $b ) = @_; my $d = eval { ( abs( $a - $b ) )->max } // 1e9; $d < 1e-12 }

# test data: N quaternions and vectors in every form
my $N = 5;
my $qp = random( 4, $N ) - 0.5; $qp /= sqrt( ( $qp**2 )->sumover )->dummy(0);
my $vp = random( 3, $N ) - 0.5;
my $tp = random($N) * 3 - 1;
my %Q = ( single => [ $qp->slice(':,(0)')->list ], array => unpdl($qp),
          columns => { map { my $i = $_; ( (qw(x y z w))[$i] => [ $qp->slice("($i)")->list ] ) } 0 .. 3 } );
my %V = ( single => [ $vp->slice(':,(0)')->list ], array => unpdl($vp),
          columns => { map { my $i = $_; ( (qw(x y z))[$i] => [ $vp->slice("($i)")->list ] ) } 0 .. 2 } );
my %T = ( single => $tp->at(0), array => [ $tp->list ] );
my $deep = sub { my $x = shift; ref $x eq 'HASH' ? { map { $_ => [ @{ $x->{$_} } ] } keys %$x }
                                : ref $x eq 'ARRAY' ? [ map { ref $_ ? [@$_] : $_ } @$x ] : $x };
my $pdl_of = sub { my ( $form, $p ) = @_; $form eq 'single' ? $p->slice(':,(0)') : $p };   # matching ndarray

# --- rotate_vec: every form of q with every form of v, new and into
for my $fq (qw( single array columns )) {
    for my $fv (qw( single array columns )) {
        my $want = rotq( $pdl_of->( $fq, $qp ), $pdl_of->( $fv, $vp ) );
        my $r = qd_rotate_vec( $Q{$fq}, $V{$fv} );
        ok( same( to_pdl( $r, [qw(x y z)] ), $want ), "rotate_vec q $fq, v $fv: new" );
        my $res_form = $fv ne 'single' ? $fv : $fq;
        my $o = $res_form eq 'columns' ? { x => [ (9) x 9 ], y => [], z => [], extra => 1 }
              : $res_form eq 'array'   ? [ map { [ 1, 2, 3, 4 ] } 1 .. 9 ] : [ 9, 9, 9, 9, 9 ];
        my $r2 = qd_rotate_vec( $Q{$fq}, $V{$fv}, $o );
        ok( refaddr($r2) == refaddr($o) && same( to_pdl( $o, [qw(x y z)] ), $want ), "rotate_vec q $fq, v $fv: into, truncated" );
        ok( $res_form ne 'array' || !grep( { @$_ != 3 } @$o ), "rotate_vec q $fq, v $fv: inner arrays truncated" );
        ok( $res_form ne 'columns' || $o->{extra} == 1, "rotate_vec q $fq, v $fv: other keys kept" );
    }
    for my $fv (qw( single array columns )) {                              # in place on v
        next if $fv eq 'single' && $fq ne 'single';                        # result is multi-valued
        my $v = $deep->( $V{$fv} );
        qd_rotate_vec( $Q{$fq}, $v, $v );
        ok( same( to_pdl( $v, [qw(x y z)] ), rotq( $pdl_of->( $fq, $qp ), $pdl_of->( $fv, $vp ) ) ), "rotate_vec q $fq, v $fv: in place" );
    }
}
ok( !eval { my $q = $deep->( $Q{array} ); qd_rotate_vec( $q, $V{array}, $q ); 1 } && $@ =~ /q cannot also be the output/, 'rotate_vec: q as output croaks' );

# --- mult and pow
for my $f1 (qw( single array columns )) {
    for my $f2 (qw( single array columns )) {
        my $want = qmult( $pdl_of->( $f1, $qp ), $pdl_of->( $f2, $qp->slice(':,-1:0') ) );
        my %Q2 = ( single => [ $qp->slice(':,(-1)')->list ], array => unpdl( $qp->slice(':,-1:0') ),
                   columns => { map { my $i = $_; ( (qw(x y z w))[$i] => [ $qp->slice("($i),-1:0")->list ] ) } 0 .. 3 } );
        ok( same( to_pdl( qd_mult( $Q{$f1}, $Q2{$f2} ), [qw(x y z w)] ), $want ), "mult $f1 x $f2: new" );
        my $a = $deep->( $Q{$f1} ); my $b = $deep->( $Q2{$f2} );
        if ( $f1 ne 'single' || $f2 eq 'single' ) {
            qd_mult( $a, $Q2{$f2}, $a ); ok( same( to_pdl( $a, [qw(x y z w)] ), $want ), "mult $f1 x $f2: in place on q1" );
        }
        if ( $f2 ne 'single' || $f1 eq 'single' ) {
            next if $f1 ne 'single' && $f1 ne $f2;                         # result takes q1's form
            qd_mult( $Q{$f1}, $b, $b ); ok( same( to_pdl( $b, [qw(x y z w)] ), $want ), "mult $f1 x $f2: in place on q2" );
        }
    }
}
{   # composition: rotating by q1 q2 is rotating by q2, then by q1
    my ( $q1, $q2 ) = ( $Q{array}, unpdl( $qp->slice(':,-1:0') ) );
    ok( same( to_pdl( qd_rotate_vec( qd_mult( $q1, $q2 ), $V{array} ), [] ),
              to_pdl( qd_rotate_vec( $q1, qd_rotate_vec( $q2, $V{array} ) ), [] ) ), 'mult: composition' );
}
for my $fq (qw( single array columns )) {
    for my $ft (qw( single array )) {
        my $want = qpow( $pdl_of->( $fq, $qp ), $ft eq 'single' ? pdl( $T{single} ) : $tp );
        ok( same( to_pdl( qd_pow( $Q{$fq}, $T{$ft} ), [qw(x y z w)] ), $want ), "pow q $fq, t $ft: new" );
    }
}
ok( same( to_pdl( qd_pow( $Q{array}, 2 ), [] ), to_pdl( qd_mult( $Q{array}, $Q{array} ), [] ) ), 'pow: q ** 2 == q * q' );
{ my $q = $deep->( $Q{array} ); qd_pow( $q, 0.5, $q );
  ok( same( to_pdl( $q, [] ), qpow( $qp, 0.5 ) ), 'pow: in place on q' ); }
ok( !eval { my $t = [ $tp->list ]; qd_pow( $Q{array}, $t, $t ); 1 }, 'pow: t as output croaks' );

# --- to_matrix
for my $fq (qw( single array )) {
    my $want = q2mat( $pdl_of->( $fq, $qp ) );
    my $m = qd_to_matrix( $Q{$fq} );
    ok( same( pdl($m), $want ), "to_matrix q $fq: new" );
    my $o = $fq eq 'single' ? [ [9], [9, 9, 9, 9] ] : [ map { [ [9] ] } 1 .. 9 ];
    qd_to_matrix( $Q{$fq}, $o );
    ok( same( pdl($o), $want ) && ( $fq eq 'single' ? @$o == 3 && !grep( { @$_ != 3 } @$o ) : @$o == $N ), "to_matrix q $fq: into, truncated" );
}
ok( same( ( pdl( qd_to_matrix( $Q{array} ) ) x $vp->dummy(0)->xchg(0,1)->dummy(0) )->clump(2)->slice('(0)'), rotq( $qp, $vp ) )
    || same( ( pdl( qd_to_matrix( $Q{array} ) ) * $vp->dummy(1) )->sumover, rotq( $qp, $vp ) ), 'to_matrix: T . v == rotate_vec' );
ok( !eval { qd_to_matrix( $Q{columns} ); 1 } && $@ =~ /no columns form/, 'to_matrix: columns croak' );
ok( !eval { my $q = $deep->( $Q{array} ); qd_to_matrix( $q, $q ); 1 }, 'to_matrix: in place croaks' );

# --- ndarrays, through the PDL side
ok( same( qd_rotate_vec( $qp, $vp ), rotq( $qp, $vp ) ), 'ndarray: rotate_vec, both ndarrays' );
ok( same( qd_rotate_vec( $Q{single}, $vp ), rotq( $pdl_of->( 'single', $qp ), $vp ) ), 'ndarray: rotate_vec, single Perl q' );
ok( same( qd_mult( $qp, $Q{single} ), qmult( $qp, $pdl_of->( 'single', $qp ) ) ), 'ndarray: mult, single Perl q2' );
ok( same( qd_pow( $qp, 0.5 ), qpow( $qp, 0.5 ) ), 'ndarray: pow, Perl scalar t' );
ok( same( qd_to_matrix($qp), q2mat($qp) ), 'ndarray: to_matrix' );
{ my $o = zeroes( 3, $N ); my $r = qd_rotate_vec( $qp, $vp, $o ); ok( refaddr($r) == refaddr($o) && same( $o, rotq( $qp, $vp ) ), 'ndarray: into' ); }
{ my $w = $vp->copy; qd_rotate_vec( $qp, $w, $w ); ok( same( $w, rotq( $qp, $vp ) ), 'ndarray: in place' ); }
{ my $o = PDL->null; qd_mult( $qp, $qp, $o ); ok( same( $o, qmult( $qp, $qp ) ), 'ndarray: null output' ); }
ok( !eval { qd_rotate_vec( $qp, $vp, zeroes( float, 3, $N ) ); 1 } && $@ =~ /double/, 'ndarray: float output croaks' );
ok( !eval { qd_rotate_vec( $Q{array}, $vp ); 1 } && $@ =~ /cannot be mixed/, 'ndarray: multi-valued Perl data with ndarray croaks' );

# --- errors
ok( !eval { qd_rotate_vec( $Q{array}, unpdl( $vp->slice(':,0:2') ) ); 1 } && $@ =~ /different lengths/, 'error: lengths differ' );
ok( !eval { qd_rotate_vec( $Q{single}, [ [ [1] ] ] ); 1 } && $@ =~ /nested/, 'error: wrong nesting' );
ok( !eval { qd_rotate_vec( $Q{single}, { x => [1], y => [1] } ); 1 } && $@ =~ /'z'/, 'error: missing column' );
ok( !eval { qd_rotate_vec( $Q{single}, 42 ); 1 } && $@ =~ /not a recognised data form/, 'error: plain number as a vector' );
ok( !eval { qd_rotate_vec( $Q{single}, [ 1, 2 ] ); 1 } && $@ =~ /must have 3 elements/, 'error: short vector' );
ok( !eval { qd_rotate_vec( $Q{single}, $V{array}, {} ); 1 } && $@ =~ /array reference/, 'error: output of the wrong kind' );
ok( !eval { qd_rotate_vec( $Q{single} ); 1 } && $@ =~ /expected 2 or 3 arguments/, 'error: too few arguments' );
ok( same( to_pdl( qd_rotate_vec( $Q{single}, [] ), [] ), zeroes( 3, 0 ) ) || !@{ qd_rotate_vec( $Q{single}, [] ) }, 'empty array of values gives an empty result' );

printf "%d/%d checks passed\n", $n - $fail, $n;
exit( $fail ? 1 : 0 );
