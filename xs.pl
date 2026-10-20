use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QRotPerl; use Time::HiRes qw(time);

# XS functions on Perl data (through the shared C code in qrot.h) vs. the
# pure-Perl implementations, plus a check that the PP kernel built on the
# same shared code is as fast as the hand-written one.

my @q = (0.1, -0.3, 0.5, 0.8);
{ my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

# correctness: XS against pure Perl, and rotq_shared against rotq
{
    my $vp = random(3, 1000) - 0.5;
    my $aoa = unpdl($vp);
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    my $maxdiff = sub { my ($a, $b) = @_; (abs(pdl($a) - pdl($b)))->max };
    printf "max |diff|: xs_rot1 %.1e  xs_rot1_list %.1e  xs_rot_aoa %.1e  xs_rot_cols %.1e  rotq_shared %.1e\n",
        $maxdiff->([ map { xs_rot1(\@q, $_) } @$aoa ], [ map { rot1(\@q, $_) } @$aoa ]),
        $maxdiff->([ map { [ xs_rot1_list(@q, @$_) ] } @$aoa ], [ map { [ rot1_list(@q, @$_) ] } @$aoa ]),
        $maxdiff->(xs_rot_aoa(\@q, $aoa), rot_aoa(\@q, $aoa)),
        $maxdiff->([ xs_rot_cols(\@q, @cols) ], [ rot_cols(\@q, @cols) ]),
        (abs(rotq_shared($qp, $vp) - rotq($qp, $vp)))->max;
}
printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;

# (1) one vector per call; 1000 calls per timed closure to amortise the closure call
{
    my $K = 1000;
    my @v = (0.3, -0.2, 0.9);
    my $vp = pdl(@v); my $rp = zeroes(3);
    my @rows = (
        ['empty loop (baseline)',              sub { for (1..$K) { } }],
        ['Perl  rot1(\@q, \@v)      -> [..]',  sub { my $r; $r = rot1(\@q, \@v) for 1..$K }],
        ['XS    xs_rot1(\@q, \@v)   -> [..]',  sub { my $r; $r = xs_rot1(\@q, \@v) for 1..$K }],
        ['Perl  rot1_list(@q, @v)   -> list',  sub { my @r; @r = rot1_list(@q, @v) for 1..$K }],
        ['XS    xs_rot1_list(@q, @v) -> list', sub { my @r; @r = xs_rot1_list(@q, @v) for 1..$K }],
        ['PDL   rotq($q, $v)  (ndarrays)',     sub { my $r; $r = rotq($qp, $vp) for 1..$K }],
        ['PDL   rotq($q, $v, $r) (prealloc)',  sub { rotq($qp, $vp, $rp) for 1..$K }],
    );
    print "(1) one vector per call: ns per call\n";
    printf "  %-38s %8.0f\n", $_->[0], 1e9 * bench($_->[1]) / $K for @rows;
}

# (2) N vectors per call
my @Ns = (1, 3, 10, 30, 100, 1e3, 1e4, 1e5, 1e6);
print "\n(2) N vectors per call: ns per vector\n";
printf "%9s %9s %9s %9s | %9s %9s %9s | %9s\n", 'N',
    'Perl aoa', 'Perl mat', 'XS aoa', 'Perl cols', 'Perl mat', 'XS cols', 'PDL rotq';
printf "%9s %29s | %29s | %9s\n", '', '--- array of [x,y,z] ---', '--- three column arrays ---', 'ndarrays';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my $aoa = unpdl($vp);
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    printf "%9.0e %9.0f %9.0f %9.0f | %9.0f %9.0f %9.0f | %9.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rot_aoa(\@q, $aoa) },
        sub { rot_aoa_mat(\@q, $aoa) },
        sub { xs_rot_aoa(\@q, $aoa) },
        sub { my @r = rot_cols(\@q, @cols) },
        sub { my @r = rot_cols_mat(\@q, @cols) },
        sub { my @r = xs_rot_cols(\@q, @cols) },
        sub { rotq($qp, $vp) };
}

# (3) PP kernel through the shared C function vs. the hand-written PP kernel
print "\n(3) PP kernels: hand-written rotq vs. rotq_shared (qrot.h): ns per vector\n";
printf "%9s %9s %11s %12s %14s\n", 'N', 'rotq', 'rotq_shared', 'rotq prealloc', 'shared prealloc';
for my $N (1e2, 1e4, 1e6, 1e7) {
    my $vp = random(3, $N) - 0.5; my $rp = zeroes(3, $N);
    printf "%9.0e %9.1f %11.1f %12.1f %14.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rotq($qp, $vp) }, sub { rotq_shared($qp, $vp) },
        sub { rotq($qp, $vp, $rp) }, sub { rotq_shared($qp, $vp, $rp) };
}
