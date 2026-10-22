use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QDispatch; use Time::HiRes qw(time);

# Dispatch of ndarrays: rotq called directly, against the Perl dispatcher
# (QDispatch::rotate_vec), the XS dispatcher handing ndarrays to the Perl one
# (xsd_rotate_vec, call_pv), and the XS dispatcher calling the PDL side's C
# function through the pointer published in $QBench::PP_ROTATE_VEC
# (xsd2_rotate_vec).  ns per call; "overhead" columns are dispatch - direct.

my @q = (0.1, -0.3, 0.5, 0.8); { my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $qp = pdl(@q);

sub bench {
    my ($code, $reps) = @_;
    $code->();
    my ($n, $t0) = (0, time);
    do { $code->(); $n++ } while time - $t0 < 0.5;
    return 1e9 * (time - $t0) / $n / $reps;
}

my %f = ( perl => \&QDispatch::rotate_vec, xsd => \&QBench::xsd_rotate_vec, xsd2 => \&QBench::xsd2_rotate_vec );

printf "perl %vd, PDL %s\n\n", $^V, $PDL::VERSION;
printf "%-32s %8s %9s %9s %9s | %8s %8s %8s\n", 'case (ns per call)', 'direct', 'Perl', 'XS call', 'XS ptr',
    'Perl ovh', 'XS ovh', 'ptr ovh';
for my $N (1, 1000) {
    my $K = $N == 1 ? 200 : 20;
    my $v = random(3, $N) - 0.5; my $o = zeroes(3, $N);
    for my $case ([ 'into, quaternion as ndarray', $qp, $o ], [ 'into, quaternion as Perl data', \@q, $o ],
                  [ 'new, quaternion as ndarray', $qp, undef ]) {
        my ($label, $qq, $out) = @$case;
        my $direct = defined $out ? sub { rotq($qp, $v, $out) for 1..$K } : sub { my $r; $r = rotq($qp, $v) for 1..$K };
        my @t = map { bench($_, $K) } $direct,
            map { my $g = $f{$_}; defined $out ? sub { $g->($qq, $v, $out) for 1..$K } : sub { my $r; $r = $g->($qq, $v) for 1..$K } }
            qw( perl xsd xsd2 );
        printf "%-32s %8.0f %9.0f %9.0f %9.0f | %8.0f %8.0f %8.0f\n", "N=$N, $label", @t, map { $_ - $t[0] } @t[1..3];
    }
}
