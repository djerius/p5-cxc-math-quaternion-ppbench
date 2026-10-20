use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QRotPerl; use Time::HiRes qw(time);

# Writing results into arrays the caller already has (the _into variants)
# vs. returning new arrays, in XS and in pure Perl.  "into" reuses a
# separate output that already holds the previous call's results; "in place"
# passes the input as the output.

my @q = (0.1, -0.3, 0.5, 0.8);
{ my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 0.5; (time-$t0)/$n }

# correctness: each _into variant against the version that returns new arrays
{
    my $vp = random(3, 1000) - 0.5;
    my $aoa = unpdl($vp);
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    my $maxdiff = sub { (abs(pdl($_[0]) - pdl($_[1])))->max };
    my (@r1, @r2, @aoa_x, @aoa_p, @cx, @cp);
    xs_rot1_into(\@q, $aoa->[0], \@r1);  rot1_into(\@q, $aoa->[0], \@r2);
    xs_rot_aoa_into(\@q, $aoa, \@aoa_x); rot_aoa_mat_into(\@q, $aoa, \@aoa_p);
    @cx = map { [] } 1..3; @cp = map { [] } 1..3;
    xs_rot_cols_into(\@q, @cols, @cx);  rot_cols_mat_into(\@q, @cols, @cp);
    my $inplace = [ map { [@$_] } @$aoa ]; xs_rot_aoa_into(\@q, $inplace, $inplace);
    printf "max |diff| vs new-array versions: xs_rot1_into %.1e  rot1_into %.1e  xs_rot_aoa_into %.1e  "
         . "rot_aoa_mat_into %.1e  xs_rot_cols_into %.1e  rot_cols_mat_into %.1e  in place %.1e\n",
        $maxdiff->(\@r1, xs_rot1(\@q, $aoa->[0])), $maxdiff->(\@r2, rot1(\@q, $aoa->[0])),
        $maxdiff->(\@aoa_x, xs_rot_aoa(\@q, $aoa)), $maxdiff->(\@aoa_p, rot_aoa_mat(\@q, $aoa)),
        $maxdiff->(\@cx, [ xs_rot_cols(\@q, @cols) ]), $maxdiff->(\@cp, [ rot_cols_mat(\@q, @cols) ]),
        $maxdiff->($inplace, xs_rot_aoa(\@q, $aoa));
}
printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;

# (1) one vector per call; 1000 calls per timed closure
{
    my $K = 1000;
    my @v = (0.3, -0.2, 0.9); my @r = (0, 0, 0);
    my @rows = (
        ['empty loop (baseline)',                 sub { for (1..$K) { } }],
        ['Perl  rot1(\@q, \@v)          -> new',   sub { my $r; $r = rot1(\@q, \@v) for 1..$K }],
        ['Perl  rot1_into(\@q, \@v, \@r)  -> into', sub { rot1_into(\@q, \@v, \@r) for 1..$K }],
        ['XS    xs_rot1(\@q, \@v)       -> new',   sub { my $r; $r = xs_rot1(\@q, \@v) for 1..$K }],
        ['XS    xs_rot1_into(\@q,\@v,\@r) -> into', sub { xs_rot1_into(\@q, \@v, \@r) for 1..$K }],
        ['XS    xs_rot1_list(@q, @v)    -> list',  sub { my @r; @r = xs_rot1_list(@q, @v) for 1..$K }],
    );
    print "(1) one vector per call: ns per call\n";
    printf "  %-40s %8.0f\n", $_->[0], 1e9 * bench($_->[1]) / $K for @rows;
}

my @Ns = (1, 10, 100, 1e3, 1e4, 1e5, 1e6);

print "\n(2) array of [x,y,z]: ns per vector\n";
printf "%9s %9s %9s | %9s %9s %9s\n", 'N', 'Perl new', 'Perl into', 'XS new', 'XS into', 'XS inplace';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my $aoa = unpdl($vp);
    my $work = [ map { [@$_] } @$aoa ];       # rotated in place, repeatedly
    my (@out_p, @out_x);
    printf "%9.0e %9.0f %9.0f | %9.0f %9.0f %9.0f\n", $N, map { 1e9 * bench($_) / $N }
        sub { rot_aoa_mat(\@q, $aoa) },
        sub { rot_aoa_mat_into(\@q, $aoa, \@out_p) },
        sub { xs_rot_aoa(\@q, $aoa) },
        sub { xs_rot_aoa_into(\@q, $aoa, \@out_x) },
        sub { xs_rot_aoa_into(\@q, $work, $work) };
}

print "\n(3) three column arrays: ns per vector\n";
printf "%9s %9s %9s | %9s %9s %9s | %12s\n", 'N', 'Perl new', 'Perl into', 'XS new', 'XS into', 'XS inplace', 'PDL prealloc';
for my $N (@Ns) {
    my $vp = random(3, $N) - 0.5;
    my @cols = map { [ $vp->slice("($_)")->list ] } 0..2;
    my @work = map { [@$_] } @cols;
    my @out_p = map { [] } 1..3; my @out_x = map { [] } 1..3;
    my @cp = map { $vp->slice("($_)")->copy } 0..2; my @op = map { zeroes($N) } 1..3;
    printf "%9.0e %9.0f %9.0f | %9.0f %9.0f %9.0f | %12.1f\n", $N, map { 1e9 * bench($_) / $N }
        sub { my @r = rot_cols_mat(\@q, @cols) },
        sub { rot_cols_mat_into(\@q, @cols, @out_p) },
        sub { my @r = xs_rot_cols(\@q, @cols) },
        sub { xs_rot_cols_into(\@q, @cols, @out_x) },
        sub { xs_rot_cols_into(\@q, @work, @work) },
        sub { rotq3($qp, @cp, @op) };
}
