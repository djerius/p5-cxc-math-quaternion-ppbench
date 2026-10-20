use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use QBench; use Time::HiRes qw(time);

# Does the memory layout of Perl scalars, left by earlier allocations, set the
# speed of XS on large Perl arrays?
#
#   perl svlayout.pl clean   # arrays built in a fresh process
#   perl svlayout.pl churn   # first create and free 3 x 3e6 scalars, as a
#                            # benchmark of the new-array API does
#
# Run each mode in its own process.

my $mode = shift // 'clean';
die "usage: $0 clean|churn\n" unless $mode =~ /^(clean|churn)$/;
my $N = 1e6;

my @q = (0.1, -0.3, 0.5, 0.8);
{ my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }

sub bench { my ($c)=@_; $c->(); my ($n,$t0)=(0,time); do { $c->(); $n++ } while time-$t0 < 1; (time-$t0)/$n }

if ($mode eq 'churn') {
    for (1..3) { my @tmp = map { [ map { rand } 1..$N ] } 1..3; my @r = xs_rot_cols(\@q, @tmp) }
}
my @in  = map { [ map { rand() - 0.5 } 1..$N ] } 1..3;
my @out = map { [] } 1..3;
xs_rot_cols_into(\@q, @in, @out);           # create the output scalars once

printf "%-6s N=%.0e  ns/vector:  into %5.1f  in place %5.1f  new %5.1f\n", $mode, $N,
    map { 1e9 * bench($_) / $N }
    sub { xs_rot_cols_into(\@q, @in, @out) },
    sub { xs_rot_cols_into(\@q, @in, @in) },
    sub { my @r = xs_rot_cols(\@q, @in) };
