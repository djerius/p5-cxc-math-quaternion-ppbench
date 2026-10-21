use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QDispatch; use Time::HiRes qw(time);

# Cost of the API's argument dispatch: one entry point (rotate_vec) that works
# out the data form of its arguments and the output mode, prototyped in pure
# Perl (QDispatch::rotate_vec) and in XS (xsd_rotate_vec), against calling the
# matching kernel directly.  ns per call; "overhead" columns are dispatch - direct.

my @q = (0.1, -0.3, 0.5, 0.8); { my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

# time $code, which makes $reps calls; ns per call
sub bench {
    my ($code, $reps) = @_;
    $code->();
    my ($n, $t0) = (0, time);
    do { $code->(); $n++ } while time - $t0 < 0.5;
    return 1e9 * (time - $t0) / $n / $reps;
}

my $PD = \&QDispatch::rotate_vec;
my $XD = \&QBench::xsd_rotate_vec;

printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;
printf "%-34s %9s %10s %10s | %9s %9s\n", 'case (ns per call)', 'direct', 'Perl disp', 'XS disp', 'Perl ovh', 'XS ovh';

sub row {
    my ($label, $reps, $direct, $perl, $xs) = @_;
    my @t = map { bench($_, $reps) } $direct, $perl, $xs;
    printf "%-34s %9.0f %10.0f %10.0f | %9.0f %9.0f\n", $label, @t, $t[1] - $t[0], $t[2] - $t[0];
}

# single vector
{
    my $K = 1000; my @v = (0.3, -0.2, 0.9); my @o = (0, 0, 0); my @w = @v;
    row('single, new',      $K, sub { my $r; $r = xs_rot1(\@q, \@v)           for 1..$K },
                                sub { my $r; $r = $PD->(\@q, \@v)              for 1..$K },
                                sub { my $r; $r = $XD->(\@q, \@v)              for 1..$K });
    row('single, into',     $K, sub { xs_rot1_into(\@q, \@v, \@o)              for 1..$K },
                                sub { $PD->(\@q, \@v, \@o)                     for 1..$K },
                                sub { $XD->(\@q, \@v, \@o)                     for 1..$K });
    row('single, in place', $K, sub { xs_rot1_into(\@q, \@w, \@w)              for 1..$K },
                                sub { $PD->(\@q, \@w, \@w)                     for 1..$K },
                                sub { $XD->(\@q, \@w, \@w)                     for 1..$K });
}

# records and columns
for my $N (10, 1000) {
    my $K = $N == 10 ? 100 : 1;
    my $vp = random(3, $N) - 0.5; my $aoa = unpdl($vp); my @ro; my $w = [ map { [@$_] } @$aoa ];
    my %col = map { my $i = $_; ( (qw(x y z))[$i] => [ $vp->slice("($i)")->list ] ) } 0..2;
    my @cin = @col{qw(x y z)}; my %co; my @cout = map { $co{$_} = [] } qw(x y z);
    row("records N=$N, new",  $K, sub { my $r; $r = xs_rot_aoa(\@q, $aoa)      for 1..$K },
                                  sub { my $r; $r = $PD->(\@q, $aoa)            for 1..$K },
                                  sub { my $r; $r = $XD->(\@q, $aoa)            for 1..$K });
    row("records N=$N, into", $K, sub { xs_rot_aoa_into(\@q, $aoa, \@ro)       for 1..$K },
                                  sub { $PD->(\@q, $aoa, \@ro)                  for 1..$K },
                                  sub { $XD->(\@q, $aoa, \@ro)                  for 1..$K });
    row("columns N=$N, into", $K, sub { xs_rot_cols_into(\@q, @cin, @cout)     for 1..$K },
                                  sub { $PD->(\@q, \%col, \%co)                 for 1..$K },
                                  sub { $XD->(\@q, \%col, \%co)                 for 1..$K });
}

# ndarrays: direct is rotq with the quaternion already an ndarray; the
# dispatchers get the quaternion as a Perl array ref and must convert it
for my $N (1, 1000) {
    my $K = $N == 1 ? 100 : 10;
    my $vp = random(3, $N) - 0.5; my $op = zeroes(3, $N);
    row("ndarray N=$N, into (q Perl)", $K, sub { rotq($qp, $vp, $op)  for 1..$K },
                                           sub { $PD->(\@q, $vp, $op) for 1..$K },
                                           sub { $XD->(\@q, $vp, $op) for 1..$K });
    row("ndarray N=$N, into (q ndarray)", $K, sub { rotq($qp, $vp, $op)  for 1..$K },
                                              sub { $PD->($qp, $vp, $op) for 1..$K },
                                              sub { $XD->($qp, $vp, $op) for 1..$K });
}
