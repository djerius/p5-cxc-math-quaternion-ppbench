use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# Pure-Perl rotation baselines vs. the PP kernels.
# Vectors in Perl are either an array of [x,y,z] refs (AoA) or three column arrays.

my @q = (0.1, -0.3, 0.5, 0.8);
{ my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

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

# correctness check against rotq
{
    my $vp = random(3, 1000) - 0.5;
    my $aoa = unpdl($vp);
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    my $ref = rotq($qp, $vp);
    my $diff = sub { (abs(pdl($_[0]) - $ref))->max };
    printf "max |diff| vs rotq: rot1 %.1e  aoa %.1e  aoa_mat %.1e  cols %.1e  cols_mat %.1e\n",
        $diff->([ map { rot1(\@q, $_) } @$aoa ]), $diff->(rot_aoa(\@q, $aoa)),
        $diff->(rot_aoa_mat(\@q, $aoa)),
        (abs(cat(map { pdl($_) } rot_cols(\@q, @cols))->xchg(0,1) - $ref))->max,
        (abs(cat(map { pdl($_) } rot_cols_mat(\@q, @cols))->xchg(0,1) - $ref))->max;
}
printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;

my @Ns = (1, 3, 10, 30, 100, 1e3, 1e4, 1e5, 1e6);

print "(1) pure Perl: ns per vector\n";
printf "%9s %9s %9s %9s %9s %9s\n", 'N', 'rot1 map', 'aoa', 'aoa mat', 'cols', 'cols mat';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my $aoa = unpdl($vp);
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    printf "%9.0e %9.0f %9.0f %9.0f %9.0f %9.0f\n", $N, map { 1e9 * bench($_) / $N }
        sub { my @r = map { rot1(\@q, $_) } @$aoa },
        sub { rot_aoa(\@q, $aoa) },
        sub { rot_aoa_mat(\@q, $aoa) },
        sub { my @r = rot_cols(\@q, @cols) },
        sub { my @r = rot_cols_mat(\@q, @cols) };
}

print "\n(2) PDL, data already in ndarrays vs. round trip from/to Perl arrays: ns per vector\n";
printf "%9s %9s %9s | %11s %11s\n", 'N', 'rotq', 'rotq3', 'AoA->rotq', 'cols->rotq3';
printf "%9s %19s | %23s\n", '', '-- PDL data --', '-- Perl data, round trip --';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my @cp = map { $vp->slice("($_)")->copy } 0..2;
    my $aoa = unpdl($vp);
    my @cols = map { [ $_->list ] } @cp;
    printf "%9.0e %9.0f %9.0f | %11.0f %11.0f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($qp, $vp) },
        sub { my @r = rotq3($qp, @cp) },
        sub { my $r = unpdl(rotq($qp, pdl($aoa))) },
        sub { my @r = map { [ $_->list ] } rotq3($qp, map { pdl($_) } @cols) };
}

print "\n(3) converting between Perl arrays and ndarrays: ns per vector\n";
printf "%9s %11s %11s | %11s %11s\n", 'N', 'pdl(AoA)', 'unpdl', 'pdl(cols)', 'cols list';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my @cp = map { $vp->slice("($_)")->copy } 0..2;
    my $aoa = unpdl($vp);
    my @cols = map { [ $_->list ] } @cp;
    printf "%9.0e %11.0f %11.0f | %11.0f %11.0f\n", $N, map { 1e9 * bench($_) / $N }
        sub { pdl($aoa) },
        sub { unpdl($vp) },
        sub { my @p = map { pdl($_) } @cols },
        sub { my @r = map { [ $_->list ] } @cp };
}
