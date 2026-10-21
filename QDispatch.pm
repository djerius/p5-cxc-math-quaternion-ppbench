package QDispatch;

# Prototype of the API's argument dispatch, in pure Perl (see dispatch.pl).
# rotate_vec( $q, $v, [$out] ) works out the data forms of $q, $v and $out,
# checks shapes, detects in place, and calls the existing kernels: the XS
# functions for Perl data and the PP kernel rotq for ndarrays.
# Prototype limits: Perl data needs a single quaternion [ $x, $y, $z, $w ].

use strict; use warnings;
use Carp qw( croak );
use Scalar::Util qw( blessed reftype refaddr );
use QBench;

use constant { SINGLE => 1, RECORDS => 2, COLUMNS => 3, NDARRAY => 4, OBJECT => 5 };

# data form of an argument, looking only at the top level and first element
sub form {
    my ( $x ) = @_;
    my $rt = reftype( $x ) // croak( 'not a reference' );
    if ( blessed $x ) {
        return NDARRAY if $x->isa( 'PDL' );
        return OBJECT;
    }
    if ( $rt eq 'ARRAY' ) {
        return RECORDS unless @$x;                      # empty batch
        return ref $x->[0] ? RECORDS : SINGLE;
    }
    return COLUMNS if $rt eq 'HASH';
    croak( 'not a recognised data form' );
}

sub rotate_vec {
    my ( $q, $v, $out ) = @_;
    my $fv = form( $v );
    my $fq = form( $q );

    if ( $fv == NDARRAY ) {
        my $qp = $fq == NDARRAY ? $q : PDL->pdl( $q );   # single quaternion as Perl data
        return QBench::rotq( $qp, $v ) unless defined $out;
        croak( 'output must be an ndarray' ) unless blessed $out && $out->isa( 'PDL' );
        croak( 'output must be double' ) unless $out->isnull || $out->type == PDL::double();
        QBench::rotq( $qp, $v, $out );                   # in place when $out is $v
        return $out;
    }

    croak( 'prototype: q must be a single [x, y, z, w]' ) unless $fq == SINGLE;

    if ( $fv == SINGLE ) {
        return QBench::xs_rot1( $q, $v ) unless defined $out;
        croak( 'output must be an array reference' ) unless ( reftype( $out ) // '' ) eq 'ARRAY';
        QBench::xs_rot1_into( $q, $v, $out );
        $#$out = 2;
        return $out;
    }

    if ( $fv == RECORDS ) {
        return QBench::xs_rot_aoa( $q, $v ) unless defined $out;
        croak( 'output must be an array reference' ) unless ( reftype( $out ) // '' ) eq 'ARRAY';
        QBench::xs_rot_aoa_into( $q, $v, $out );         # in place when $out is $v
        $#$out = $#$v;     # inner records: truncating them belongs in the kernel loop
        return $out;
    }

    # columns
    my @in = @{$v}{qw( x y z )};
    croak( 'columns need x, y and z array references' )
      if grep { ( reftype( $_ ) // '' ) ne 'ARRAY' } @in;
    unless ( defined $out ) {
        my @r = QBench::xs_rot_cols( $q, @in );
        return { x => $r[0], y => $r[1], z => $r[2] };
    }
    croak( 'output must be a hash reference' ) unless ( reftype( $out ) // '' ) eq 'HASH';
    my @o = map { $out->{$_} //= [] } qw( x y z );
    QBench::xs_rot_cols_into( $q, @in, @o );         # in place when the arrays are the same
    $#$_ = $#{ $in[0] } for @o;
    return $out;
}

1;
