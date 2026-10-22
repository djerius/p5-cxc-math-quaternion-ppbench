package QDispatch;

# Prototype of the API's argument dispatch, in pure Perl (see dispatch.pl).
# rotate_vec( $q, $v, [$out] ) works out the data form of $v and the output
# mode, checks shapes, and calls the existing kernels: the XS functions for
# Perl data and the PP kernel rotq for ndarrays.
# Prototype limits: Perl data needs a single quaternion [ $x, $y, $z, $w ].
#
# Written for speed, since the pure-Perl module has to dispatch in Perl
# (see qd-variants.pl for the alternatives timed):
#   - the data form comes from one ref() per argument: 'ARRAY' or 'HASH'
#     for Perl data, '' for a non-reference, a class name for an object;
#   - it branches on that directly, most common form first, rather than
#     classifying first and branching on the result;
#   - the ndarray output check compares integer type codes, not PDL Type
#     objects.

use strict; use warnings; use feature 'state';
use Carp qw( croak );
use QBench;

sub rotate_vec {
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
        # array of values; an empty array is an empty array of values
        return QBench::xs_rot_aoa( $q, $v ) unless defined $out;
        croak( 'output must be an array reference' ) unless ref $out eq 'ARRAY';
        QBench::xs_rot_aoa_into( $q, $v, $out );     # in place when $out is $v
        $#$out = $#$v;     # inner records: truncating them belongs in the kernel loop
        return $out;
    }

    if ( $rv eq 'HASH' ) {                               # columns
        croak( 'prototype: q must be a single [x, y, z, w]' ) unless ref $q eq 'ARRAY' && !ref $q->[0];
        my ( $x, $y, $z ) = ( $v->{x}, $v->{y}, $v->{z} );
        croak( 'columns need x, y and z array references' )
          unless ref $x eq 'ARRAY' && ref $y eq 'ARRAY' && ref $z eq 'ARRAY';
        if ( !defined $out ) {
            my @r = QBench::xs_rot_cols( $q, $x, $y, $z );
            return { x => $r[0], y => $r[1], z => $r[2] };
        }
        croak( 'output must be a hash reference' ) unless ref $out eq 'HASH';
        # in place when the output arrays are the input arrays
        QBench::xs_rot_cols_into( $q, $x, $y, $z, $out->{x} //= [], $out->{y} //= [], $out->{z} //= [] );
        $#{ $out->{x} } = $#{ $out->{y} } = $#{ $out->{z} } = $#$x;
        return $out;
    }

    croak( 'not a reference' ) if $rv eq '';
    croak( 'not a recognised data form' ) unless $v->isa( 'PDL' );   # ndarray
    state $double = PDL::double()->enum;
    $q = PDL->pdl( $q ) if ref $q eq 'ARRAY';           # single quaternion as Perl data
    return QBench::rotq( $q, $v ) unless defined $out;
    croak( 'output must be a double ndarray' )
      unless ref $out && $out->isa( 'PDL' ) && ( $out->isnull || $out->get_datatype == $double );
    QBench::rotq( $q, $v, $out );                       # in place when $out is $v
    return $out;
}

1;
