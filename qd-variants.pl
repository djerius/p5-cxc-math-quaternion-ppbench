# Candidate fast versions of QDispatch::rotate_vec, timed against the
# current QDispatch.pm in one process.  ns per call, best of 5 rounds.
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QDispatch; use Carp qw( croak ); use Time::HiRes qw( time );

my $D = PDL::double()->enum;

# (a) nested blocks; one ref() per argument; branch directly on its result
sub rv_blocks {
    my ( $q, $v, $out ) = @_;
    my $rv = ref $v;
    if ( $rv eq 'ARRAY' ) {
        croak( 'prototype: q must be a single [x, y, z, w]' ) unless ref $q eq 'ARRAY' && !ref $q->[0];
        if ( @$v && !ref $v->[0] ) {                    # single value
            return QBench::xs_rot1( $q, $v ) unless defined $out;
            croak( 'output must be an array reference' ) unless ref $out eq 'ARRAY';
            QBench::xs_rot1_into( $q, $v, $out );
            $#$out = 2;
            return $out;
        }
        return QBench::xs_rot_aoa( $q, $v ) unless defined $out;   # array of values
        croak( 'output must be an array reference' ) unless ref $out eq 'ARRAY';
        QBench::xs_rot_aoa_into( $q, $v, $out );
        $#$out = $#$v;
        return $out;
    }
    if ( $rv eq 'HASH' ) {                               # columns
        croak( 'prototype: q must be a single [x, y, z, w]' ) unless ref $q eq 'ARRAY' && !ref $q->[0];
        my @in = @{$v}{qw( x y z )};
        croak( 'columns need x, y and z array references' )
          unless ref $in[0] eq 'ARRAY' && ref $in[1] eq 'ARRAY' && ref $in[2] eq 'ARRAY';
        my @r = defined $out ? () : QBench::xs_rot_cols( $q, @in );
        return { x => $r[0], y => $r[1], z => $r[2] } unless defined $out;
        croak( 'output must be a hash reference' ) unless ref $out eq 'HASH';
        my @o = ( $out->{x} //= [], $out->{y} //= [], $out->{z} //= [] );
        QBench::xs_rot_cols_into( $q, @in, @o );
        $#{ $o[0] } = $#{ $o[1] } = $#{ $o[2] } = $#{ $in[0] };
        return $out;
    }
    croak( 'not a reference' ) if $rv eq '';
    croak( 'not a recognised data form' ) unless $v->isa( 'PDL' );   # ndarray
    $q = PDL->pdl( $q ) if ref $q eq 'ARRAY';
    return QBench::rotq( $q, $v ) unless defined $out;
    croak( 'output must be a double ndarray' )
      unless ref $out && $out->isa( 'PDL' ) && ( $out->isnull || $out->get_datatype == $D );
    QBench::rotq( $q, $v, $out );
    return $out;
}

# (b) flat: no multi-statement blocks, only statement modifiers
sub rv_flat {
    my ( $q, $v, $out ) = @_;
    my $rv = ref $v;
    my $perl = $rv eq 'ARRAY' || $rv eq 'HASH';
    croak( 'prototype: q must be a single [x, y, z, w]' ) if $perl && !( ref $q eq 'ARRAY' && !ref $q->[0] );
    my $single = $rv eq 'ARRAY' && @$v && !ref $v->[0];

    return QBench::xs_rot1( $q, $v ) if $single && !defined $out;
    return QBench::xs_rot_aoa( $q, $v ) if $rv eq 'ARRAY' && !defined $out;
    croak( 'output must be an array reference' ) if $rv eq 'ARRAY' && ref $out ne 'ARRAY';
    return ( QBench::xs_rot1_into( $q, $v, $out ), $#$out = 2, $out )[-1] if $single;
    return ( QBench::xs_rot_aoa_into( $q, $v, $out ), $#$out = $#$v, $out )[-1] if $rv eq 'ARRAY';

    my @in = $rv eq 'HASH' ? @{$v}{qw( x y z )} : ();
    croak( 'columns need x, y and z array references' )
      if $rv eq 'HASH' && !( ref $in[0] eq 'ARRAY' && ref $in[1] eq 'ARRAY' && ref $in[2] eq 'ARRAY' );
    my @r = $rv eq 'HASH' && !defined $out ? QBench::xs_rot_cols( $q, @in ) : ();
    return { x => $r[0], y => $r[1], z => $r[2] } if $rv eq 'HASH' && !defined $out;
    croak( 'output must be a hash reference' ) if $rv eq 'HASH' && ref $out ne 'HASH';
    my @o = $rv eq 'HASH' ? ( $out->{x} //= [], $out->{y} //= [], $out->{z} //= [] ) : ();
    return ( QBench::xs_rot_cols_into( $q, @in, @o ), $#{ $o[0] } = $#{ $o[1] } = $#{ $o[2] } = $#{ $in[0] }, $out )[-1]
      if $rv eq 'HASH';

    croak( 'not a reference' ) if $rv eq '';
    croak( 'not a recognised data form' ) unless $v->isa( 'PDL' );
    $q = PDL->pdl( $q ) if ref $q eq 'ARRAY';
    return QBench::rotq( $q, $v ) unless defined $out;
    croak( 'output must be a double ndarray' )
      unless ref $out && $out->isa( 'PDL' ) && ( $out->isnull || $out->get_datatype == $D );
    return ( QBench::rotq( $q, $v, $out ), $out )[-1];
}

# (c) as (a), without copying @_, and with single hash fetches for columns
sub rv_lean {
    my $rv = ref $_[1];
    if ( $rv eq 'ARRAY' ) {
        croak( 'prototype: q must be a single [x, y, z, w]' ) unless ref $_[0] eq 'ARRAY' && !ref $_[0][0];
        if ( @{ $_[1] } && !ref $_[1][0] ) {                # single value
            return QBench::xs_rot1( $_[0], $_[1] ) unless defined $_[2];
            croak( 'output must be an array reference' ) unless ref $_[2] eq 'ARRAY';
            QBench::xs_rot1_into( $_[0], $_[1], $_[2] );
            $#{ $_[2] } = 2;
            return $_[2];
        }
        return QBench::xs_rot_aoa( $_[0], $_[1] ) unless defined $_[2];   # array of values
        croak( 'output must be an array reference' ) unless ref $_[2] eq 'ARRAY';
        QBench::xs_rot_aoa_into( $_[0], $_[1], $_[2] );
        $#{ $_[2] } = $#{ $_[1] };
        return $_[2];
    }
    if ( $rv eq 'HASH' ) {                               # columns
        croak( 'prototype: q must be a single [x, y, z, w]' ) unless ref $_[0] eq 'ARRAY' && !ref $_[0][0];
        my ( $x, $y, $z ) = ( $_[1]{x}, $_[1]{y}, $_[1]{z} );
        croak( 'columns need x, y and z array references' )
          unless ref $x eq 'ARRAY' && ref $y eq 'ARRAY' && ref $z eq 'ARRAY';
        if ( !defined $_[2] ) {
            my @r = QBench::xs_rot_cols( $_[0], $x, $y, $z );
            return { x => $r[0], y => $r[1], z => $r[2] };
        }
        croak( 'output must be a hash reference' ) unless ref $_[2] eq 'HASH';
        my $o = $_[2];
        QBench::xs_rot_cols_into( $_[0], $x, $y, $z, $o->{x} //= [], $o->{y} //= [], $o->{z} //= [] );
        $#{ $o->{x} } = $#{ $o->{y} } = $#{ $o->{z} } = $#$x;
        return $o;
    }
    croak( 'not a reference' ) if $rv eq '';
    croak( 'not a recognised data form' ) unless $_[1]->isa( 'PDL' );   # ndarray
    my $q = ref $_[0] eq 'ARRAY' ? PDL->pdl( $_[0] ) : $_[0];
    return QBench::rotq( $q, $_[1] ) unless defined $_[2];
    croak( 'output must be a double ndarray' )
      unless ref $_[2] && $_[2]->isa( 'PDL' ) && ( $_[2]->isnull || $_[2]->get_datatype == $D );
    QBench::rotq( $q, $_[1], $_[2] );
    return $_[2];
}

# --- check: every variant agrees with the current QDispatch on every case
my @q = ( 0.1, -0.3, 0.5, 0.8 ); { my $n = sqrt( eval { my $s = 0; $s += $_**2 for @q; $s } ); $_ /= $n for @q }
my $qp = pdl(@q);
my $vp = random( 3, 10 ) - 0.5;
my $aoa = unpdl($vp);
my %col = ( x => [ $vp->slice('(0)')->list ], y => [ $vp->slice('(1)')->list ], z => [ $vp->slice('(2)')->list ] );
my %f = ( current => \&QDispatch::rotate_vec, blocks => \&rv_blocks, flat => \&rv_flat, lean => \&rv_lean );
my $flat = sub { my $r = shift; ref $r eq 'HASH' ? pdl( @{$r}{qw(x y z)} ) : pdl($r) };
for my $name (qw( blocks flat lean )) {
    my $f = $f{$name};
    my @checks = (
        [ 'single new',        sub { $_[0]->( \@q, $aoa->[0] ) } ],
        [ 'single into',       sub { my @o = (9) x 5; $_[0]->( \@q, $aoa->[0], \@o ) } ],
        [ 'values new',        sub { $_[0]->( \@q, $aoa ) } ],
        [ 'values into',       sub { my @o; $_[0]->( \@q, $aoa, \@o ) } ],
        [ 'empty values',      sub { $_[0]->( \@q, [] ) } ],
        [ 'columns new',       sub { $_[0]->( \@q, \%col ) } ],
        [ 'columns into',      sub { my %o; $_[0]->( \@q, \%col, \%o ) } ],
        [ 'ndarray new',       sub { $_[0]->( \@q, $vp ) } ],
        [ 'ndarray into',      sub { my $o = zeroes( 3, 10 ); $_[0]->( $qp, $vp, $o ) } ],
    );
    for my $c (@checks) {
        my ( $want, $got ) = map { $c->[1]->($_) } $f{current}, $f;
        my $d = ( ref $want eq 'ARRAY' && !@$want ) ? ( @$got ? 1 : 0 ) : ( abs( $flat->($want) - $flat->($got) ) )->max;
        die "$name disagrees with current on $c->[0] ($d)\n" if $d > 1e-12;
    }
    for my $bad ( [ \@q, \"s" ], [ \@q, 42 ], [ \@q, $aoa, {} ], [ \@q, \%col, [] ], [ \@q, $vp, zeroes( float, 3, 10 ) ] ) {
        die "$name does not croak on a bad call\n" if eval { $f->(@$bad); 1 };
    }
}
print "blocks, flat and lean agree with the current version on all cases, and croak on bad calls\n\n";

# --- time
sub bench { my ( $c, $k ) = @_; $c->(); my ( $n, $t0 ) = ( 0, time ); do { $c->(); $n++ } while time - $t0 < 0.5; 1e9 * ( time - $t0 ) / $n / $k }
my $K = 2000;
my @v1 = ( 0.3, -0.2, 0.9 ); my @o1 = ( 0, 0, 0 ); my @ro; my %co; my $v1p = pdl( 0.3, -0.2, 0.9 ); my $o1p = zeroes(3);
my @cases = (
    [ 'single value, new',  sub { my $f = shift; my $r; $r = $f->( \@q, \@v1 )        for 1 .. $K } ],
    [ 'single value, into', sub { my $f = shift; $f->( \@q, \@v1, \@o1 )              for 1 .. $K } ],
    [ '10 values, into',    sub { my $f = shift; $f->( \@q, $aoa, \@ro )              for 1 .. $K } ],
    [ '10 columns, new',    sub { my $f = shift; my $r; $r = $f->( \@q, \%col )       for 1 .. $K } ],
    [ '10 columns, into',   sub { my $f = shift; $f->( \@q, \%col, \%co )             for 1 .. $K } ],
    [ 'ndarray, 1, into',   sub { my $f = shift; $f->( $qp, $v1p, $o1p )              for 1 .. $K } ],
);
my %best;
for ( 1 .. 5 ) {
    for my $c (@cases) {
        for my $name (qw( current blocks flat lean )) {
            my $t = bench( sub { $c->[1]->( $f{$name} ) }, $K );
            $best{ $c->[0] }{$name} = $t if !defined $best{ $c->[0] }{$name} || $t < $best{ $c->[0] }{$name};
        }
    }
}
printf "%-20s %10s %10s %10s %10s\n", 'case (ns per call)', 'current', 'blocks', 'flat', 'lean';
printf "%-20s %7.0f ns %7.0f ns %7.0f ns %7.0f ns\n", $_->[0], @{ $best{ $_->[0] } }{qw( current blocks flat lean )} for @cases;
