use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use Scalar::Util qw(blessed reftype refaddr); use Time::HiRes qw(time);

# Cost of the individual steps of the pure-Perl dispatcher, ns per operation.
my @q = (0.1, -0.3, 0.5, 0.8); my $qp = pdl(@q); my $v = [ 0.3, -0.2, 0.9 ]; my $aoa = [ $v, $v ];
my %col = ( x => [1], y => [2], z => [3] ); my $op = zeroes(3);
my $D = PDL::double()->enum;
sub bench { my ($c, $k) = @_; $c->(); my ($n, $t0) = (0, time); do { $c->(); $n++ } while time - $t0 < 0.5; 1e9 * (time - $t0) / $n / $k }
my $K = 10000;
my @rows = (
  [ 'empty loop',                               sub { for (1..$K) { } } ],
  [ 'call an empty Perl sub with 3 args',        sub { my $f = sub { my ($a, $b, $c) = @_ }; $f->(\@q, $v, $v) for 1..$K } ],
  [ 'reftype + blessed on an array ref',         sub { for (1..$K) { my $r = reftype($v); my $b = blessed $v } } ],
  [ 'records/single test (ref $x->[0])',         sub { for (1..$K) { my $s = ref $aoa->[0] ? 1 : 0 } } ],
  [ 'refaddr($out) == refaddr($in)',             sub { for (1..$K) { my $s = refaddr($v) == refaddr($aoa) } } ],
  [ 'hash slice @{$h}{qw(x y z)}',               sub { for (1..$K) { my @c = @col{qw(x y z)} } } ],
  [ 'blessed + ->isa("PDL") on an ndarray',      sub { for (1..$K) { my $b = blessed($qp) && $qp->isa('PDL') } } ],
  [ '$out->isnull',                              sub { for (1..$K) { my $b = $op->isnull } } ],
  [ '$out->type == PDL::double()',               sub { for (1..$K) { my $b = $op->type == PDL::double() } } ],
  [ '$out->get_datatype == $D (integer)',        sub { for (1..$K) { my $b = $op->get_datatype == $D } } ],
  [ 'PDL->pdl([4 numbers])',                     sub { for (1..$K/10) { my $p = PDL->pdl(\@q) } } ],
);
printf "%-44s %8s\n", 'step', 'ns';
for (@rows) { my $k = $_->[0] =~ /PDL->pdl/ ? $K/10 : $K; printf "%-44s %8.0f\n", $_->[0], bench($_->[1], $k) }
