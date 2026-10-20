package QRotPerl;

# Pure-Perl rotation of 3-vectors by a scalar-last unit quaternion [x, y, z, w]
# (Hamilton convention), used by pureperl.pl and xs.pl.

use strict; use warnings;
use Exporter 'import';
our @EXPORT = qw( rot1 rot1_list rot_aoa matrix rot_aoa_mat rot_cols rot_cols_mat
                  rot1_into rot_aoa_mat_into rot_cols_mat_into );

# one vector per call: the scalar API
sub rot1 {
    my ($q, $v) = @_;
    my ($x, $y, $z, $w) = @$q;
    my ($vx, $vy, $vz) = @$v;
    my $tx = 2*($y*$vz - $z*$vy); my $ty = 2*($z*$vx - $x*$vz); my $tz = 2*($x*$vy - $y*$vx);
    return [ $vx + $w*$tx + ($y*$tz - $z*$ty),
             $vy + $w*$ty + ($z*$tx - $x*$tz),
             $vz + $w*$tz + ($x*$ty - $y*$tx) ];
}

# AoA, quaternion form, loop inlined
sub rot_aoa {
    my ($q, $vs) = @_;
    my ($x, $y, $z, $w) = @$q;
    my @r;
    for my $v (@$vs) {
        my ($vx, $vy, $vz) = @$v;
        my $tx = 2*($y*$vz - $z*$vy); my $ty = 2*($z*$vx - $x*$vz); my $tz = 2*($x*$vy - $y*$vx);
        push @r, [ $vx + $w*$tx + ($y*$tz - $z*$ty),
                   $vy + $w*$ty + ($z*$tx - $x*$tz),
                   $vz + $w*$tz + ($x*$ty - $y*$tx) ];
    }
    return \@r;
}

sub matrix {
    my ($x, $y, $z, $w) = @{$_[0]};
    return ( 1-2*($y*$y+$z*$z), 2*($x*$y-$z*$w),   2*($x*$z+$y*$w),
             2*($x*$y+$z*$w),   1-2*($x*$x+$z*$z), 2*($y*$z-$x*$w),
             2*($x*$z-$y*$w),   2*($y*$z+$x*$w),   1-2*($x*$x+$y*$y) );
}

# AoA, rotation matrix computed once per call
sub rot_aoa_mat {
    my ($q, $vs) = @_;
    my ($m0,$m1,$m2,$m3,$m4,$m5,$m6,$m7,$m8) = matrix($q);
    my @r;
    for my $v (@$vs) {
        my ($vx, $vy, $vz) = @$v;
        push @r, [ $m0*$vx + $m1*$vy + $m2*$vz,
                   $m3*$vx + $m4*$vy + $m5*$vz,
                   $m6*$vx + $m7*$vy + $m8*$vz ];
    }
    return \@r;
}

# three columns in, three columns out, quaternion form
sub rot_cols {
    my ($q, $cx, $cy, $cz) = @_;
    my ($x, $y, $z, $w) = @$q;
    my (@rx, @ry, @rz);
    for my $i (0 .. $#$cx) {
        my ($vx, $vy, $vz) = ($cx->[$i], $cy->[$i], $cz->[$i]);
        my $tx = 2*($y*$vz - $z*$vy); my $ty = 2*($z*$vx - $x*$vz); my $tz = 2*($x*$vy - $y*$vx);
        push @rx, $vx + $w*$tx + ($y*$tz - $z*$ty);
        push @ry, $vy + $w*$ty + ($z*$tx - $x*$tz);
        push @rz, $vz + $w*$tz + ($x*$ty - $y*$tx);
    }
    return (\@rx, \@ry, \@rz);
}

# three columns, rotation matrix
sub rot_cols_mat {
    my ($q, $cx, $cy, $cz) = @_;
    my ($m0,$m1,$m2,$m3,$m4,$m5,$m6,$m7,$m8) = matrix($q);
    my (@rx, @ry, @rz);
    for my $i (0 .. $#$cx) {
        my ($vx, $vy, $vz) = ($cx->[$i], $cy->[$i], $cz->[$i]);
        push @rx, $m0*$vx + $m1*$vy + $m2*$vz;
        push @ry, $m3*$vx + $m4*$vy + $m5*$vz;
        push @rz, $m6*$vx + $m7*$vy + $m8*$vz;
    }
    return (\@rx, \@ry, \@rz);
}


# one vector per call, plain scalars in and out: (x, y, z, w, vx, vy, vz) -> (x', y', z')
sub rot1_list {
    my ($x, $y, $z, $w, $vx, $vy, $vz) = @_;
    my $tx = 2*($y*$vz - $z*$vy); my $ty = 2*($z*$vx - $x*$vz); my $tz = 2*($x*$vy - $y*$vx);
    return ( $vx + $w*$tx + ($y*$tz - $z*$ty),
             $vy + $w*$ty + ($z*$tx - $x*$tz),
             $vz + $w*$tz + ($x*$ty - $y*$tx) );
}


# The _into variants assign to elements of arrays the caller supplies, so
# existing scalars (and inner [x, y, z] arrays) are reused rather than
# replaced; an output may be the input itself.

sub rot1_into {
    my ($q, $v, $r) = @_;
    my ($x, $y, $z, $w) = @$q;
    my ($vx, $vy, $vz) = @$v;
    my $tx = 2*($y*$vz - $z*$vy); my $ty = 2*($z*$vx - $x*$vz); my $tz = 2*($x*$vy - $y*$vx);
    $r->[0] = $vx + $w*$tx + ($y*$tz - $z*$ty);
    $r->[1] = $vy + $w*$ty + ($z*$tx - $x*$tz);
    $r->[2] = $vz + $w*$tz + ($x*$ty - $y*$tx);
    return;
}

sub rot_aoa_mat_into {
    my ($q, $vs, $out) = @_;
    my ($m0,$m1,$m2,$m3,$m4,$m5,$m6,$m7,$m8) = matrix($q);
    for my $i (0 .. $#$vs) {
        my ($vx, $vy, $vz) = @{ $vs->[$i] };
        my $r = $out->[$i] //= [];
        $r->[0] = $m0*$vx + $m1*$vy + $m2*$vz;
        $r->[1] = $m3*$vx + $m4*$vy + $m5*$vz;
        $r->[2] = $m6*$vx + $m7*$vy + $m8*$vz;
    }
    return;
}

sub rot_cols_mat_into {
    my ($q, $cx, $cy, $cz, $ox, $oy, $oz) = @_;
    my ($m0,$m1,$m2,$m3,$m4,$m5,$m6,$m7,$m8) = matrix($q);
    for my $i (0 .. $#$cx) {
        my ($vx, $vy, $vz) = ($cx->[$i], $cy->[$i], $cz->[$i]);
        $ox->[$i] = $m0*$vx + $m1*$vy + $m2*$vz;
        $oy->[$i] = $m3*$vx + $m4*$vy + $m5*$vz;
        $oz->[$i] = $m6*$vx + $m7*$vy + $m8*$vz;
    }
    return;
}

1;
