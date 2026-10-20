use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use Time::HiRes qw(time);
use QRotPerl;

# Pure-Perl rotation baselines vs. the PP kernels.
# Vectors in Perl are either an array of [x,y,z] refs (AoA) or three column arrays.

my @q = (0.1, -0.3, 0.5, 0.8);
{ my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

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
