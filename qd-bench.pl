use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use Time::HiRes qw(time);

# The generic dispatch engine (qd_*, qdispatch.h) against the hand-written
# XS dispatcher (xsd_rotate_vec for Perl data, xsd2_rotate_vec for ndarrays)
# and the kernels called directly.  ns per call.

my @q = (0.1, -0.3, 0.5, 0.8); { my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench {
    my ($code, $reps) = @_;
    $code->();
    my ($n, $t0) = (0, time);
    do { $code->(); $n++ } while time - $t0 < 0.5;
    return 1e9 * (time - $t0) / $n / $reps;
}
sub row { my ($label, $reps, @subs) = @_; printf "%-30s" . " %10s" x @subs . "\n", $label, map { defined $_ ? sprintf('%7.0f ns', bench($_, $reps)) : '' } @subs }

printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;

print "(1) rotate_vec on Perl data\n";
printf "%-30s %10s %10s %10s\n", 'case (ns per call)', 'kernel', 'xsd', 'qd';
{
    my $K = 1000; my @v = (0.3, -0.2, 0.9); my @o = (0, 0, 0); my @w = @v;
    row('single, new', $K,      sub { my $r; $r = xs_rot1(\@q, \@v) for 1..$K }, sub { my $r; $r = xsd_rotate_vec(\@q, \@v) for 1..$K },
                                sub { my $r; $r = qd_rotate_vec(\@q, \@v) for 1..$K });
    row('single, into', $K,     sub { xs_rot1_into(\@q, \@v, \@o) for 1..$K }, sub { xsd_rotate_vec(\@q, \@v, \@o) for 1..$K },
                                sub { qd_rotate_vec(\@q, \@v, \@o) for 1..$K });
    row('single, in place', $K, sub { xs_rot1_into(\@q, \@w, \@w) for 1..$K }, sub { xsd_rotate_vec(\@q, \@w, \@w) for 1..$K },
                                sub { qd_rotate_vec(\@q, \@w, \@w) for 1..$K });
}
for my $N (10, 1000) {
    my $K = $N == 10 ? 100 : 1;
    my $vp = random(3, $N) - 0.5; my $aoa = unpdl($vp); my @ro; my $w = [ map { [@$_] } @$aoa ];
    my %col = map { my $i = $_; ( (qw(x y z))[$i] => [ $vp->slice("($i)")->list ] ) } 0..2;
    my @cin = @col{qw(x y z)}; my %co; my @cout = map { $co{$_} = [] } qw(x y z);
    my %cw = map { $_ => [ @{ $col{$_} } ] } qw(x y z);
    row("$N values, new", $K,     sub { my $r; $r = xs_rot_aoa(\@q, $aoa) for 1..$K }, sub { my $r; $r = xsd_rotate_vec(\@q, $aoa) for 1..$K },
                                  sub { my $r; $r = qd_rotate_vec(\@q, $aoa) for 1..$K });
    row("$N values, into", $K,    sub { xs_rot_aoa_into(\@q, $aoa, \@ro) for 1..$K }, sub { xsd_rotate_vec(\@q, $aoa, \@ro) for 1..$K },
                                  sub { qd_rotate_vec(\@q, $aoa, \@ro) for 1..$K });
    row("$N values, in place", $K, sub { xs_rot_aoa_into(\@q, $w, $w) for 1..$K }, sub { xsd_rotate_vec(\@q, $w, $w) for 1..$K },
                                  sub { qd_rotate_vec(\@q, $w, $w) for 1..$K });
    row("$N columns, new", $K,    sub { my @r = xs_rot_cols(\@q, @cin) for 1..$K }, sub { my $r; $r = xsd_rotate_vec(\@q, \%col) for 1..$K },
                                  sub { my $r; $r = qd_rotate_vec(\@q, \%col) for 1..$K });
    row("$N columns, into", $K,   sub { xs_rot_cols_into(\@q, @cin, @cout) for 1..$K }, sub { xsd_rotate_vec(\@q, \%col, \%co) for 1..$K },
                                  sub { qd_rotate_vec(\@q, \%col, \%co) for 1..$K });
    row("$N columns, in place", $K, sub { xs_rot_cols_into(\@q, @cw{qw(x y z)}, @cw{qw(x y z)}) for 1..$K },
                                  sub { xsd_rotate_vec(\@q, \%cw, \%cw) for 1..$K }, sub { qd_rotate_vec(\@q, \%cw, \%cw) for 1..$K });
}

print "\n(2) rotate_vec on ndarrays\n";
printf "%-30s %10s %10s %10s\n", 'case (ns per call)', 'rotq', 'xsd2', 'qd';
for my $N (1, 1000) {
    my $K = $N == 1 ? 200 : 20;
    my $vp = random(3, $N) - 0.5; my $o = zeroes(3, $N);
    row("N=$N, into, q ndarray", $K,   sub { rotq($qp, $vp, $o) for 1..$K }, sub { xsd2_rotate_vec($qp, $vp, $o) for 1..$K },
                                        sub { qd_rotate_vec($qp, $vp, $o) for 1..$K });
    row("N=$N, into, q Perl", $K,      undef, sub { xsd2_rotate_vec(\@q, $vp, $o) for 1..$K }, sub { qd_rotate_vec(\@q, $vp, $o) for 1..$K });
    row("N=$N, new, q ndarray", $K,    sub { my $r; $r = rotq($qp, $vp) for 1..$K }, sub { my $r; $r = xsd2_rotate_vec($qp, $vp) for 1..$K },
                                        sub { my $r; $r = qd_rotate_vec($qp, $vp) for 1..$K });
}

print "\n(3) the other operations through the engine\n";
printf "%-30s %10s %10s\n", 'case (ns per call)', 'PP kernel', 'qd';
{
    my $K = 1000; my @t10 = map { $_ / 10 } 1..10; my $aq = [ map { [@q] } 1..10 ];
    my $q10 = pdl($aq); my $q1000 = random(4, 1000) - 0.5; $q1000 /= sqrt(($q1000**2)->sumover)->dummy(0);
    row('mult, single x single', $K, undef, sub { my $r; $r = qd_mult(\@q, \@q) for 1..$K });
    row('mult, 10 values x single', 100, undef, sub { my $r; $r = qd_mult($aq, \@q) for 1..100 });
    row('pow, single, t number', $K, undef, sub { my $r; $r = qd_pow(\@q, 0.5) for 1..$K });
    row('pow, single, 10 t values', 100, undef, sub { my $r; $r = qd_pow(\@q, \@t10) for 1..100 });
    row('to_matrix, single', $K, undef, sub { my $r; $r = qd_to_matrix(\@q) for 1..$K });
    row('to_matrix, 10 values', 100, undef, sub { my $r; $r = qd_to_matrix($aq) for 1..100 });
    row('mult, ndarray N=1000', 20, sub { my $r; $r = qmult($q1000, $q1000) for 1..20 }, sub { my $r; $r = qd_mult($q1000, $q1000) for 1..20 });
    row('pow, ndarray N=1000', 20, sub { my $r; $r = qpow($q1000, 0.5) for 1..20 }, sub { my $r; $r = qd_pow($q1000, 0.5) for 1..20 });
    row('to_matrix, ndarray N=1000', 20, sub { my $r; $r = q2mat($q1000) for 1..20 }, sub { my $r; $r = qd_to_matrix($q1000) for 1..20 });
}
