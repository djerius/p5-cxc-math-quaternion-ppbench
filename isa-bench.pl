use v5.36;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use Scalar::Util qw( blessed reftype ); use Time::HiRes qw( time );

# Is the infix isa operator faster than reftype + blessed (+ ->isa) for the
# classification in QDispatch::form?  ns per operation, loop overhead included.
my $p   = pdl( 1, 2, 3 );
my $a   = [ 1, 2, 3 ];
my $aoa = [ [ 1, 2, 3 ] ];
my $h   = { x => [1] };
my $K   = 100_000;
sub bench { my ($c) = @_; $c->(); my ($n, $t0) = (0, time); do { $c->(); $n++ } while time - $t0 < 0.5; 1e9 * (time - $t0) / $n / $K }

# form() as in QDispatch, and with the isa operator first
sub form_now { my ($x) = @_; my $rt = reftype($x) // die; if (blessed $x) { return 4 if $x->isa('PDL'); return 5 }
               if ($rt eq 'ARRAY') { return 2 unless @$x; return ref $x->[0] ? 2 : 1 } return 3 if $rt eq 'HASH'; die }
sub form_isa { my ($x) = @_; return 4 if $x isa PDL; my $rt = reftype($x) // die; return 5 if blessed $x;
               if ($rt eq 'ARRAY') { return 2 unless @$x; return ref $x->[0] ? 2 : 1 } return 3 if $rt eq 'HASH'; die }

my @rows = (
  [ 'empty loop',                                  sub { for (1..$K) { } } ],
  [ 'ndarray: blessed($x) && $x->isa("PDL")',      sub { for (1..$K) { my $r = blessed($p) && $p->isa('PDL') } } ],
  [ 'ndarray: $x isa PDL',                         sub { for (1..$K) { my $r = $p isa PDL } } ],
  [ 'array ref: reftype($x) + blessed($x)',        sub { for (1..$K) { my $t = reftype($a); my $b = blessed $a } } ],
  [ 'array ref: blessed($x) alone',                sub { for (1..$K) { my $b = blessed $a } } ],
  [ 'array ref: $x isa PDL (false)',               sub { for (1..$K) { my $r = $a isa PDL } } ],
  [ 'array ref: ref($x) eq "ARRAY"',               sub { for (1..$K) { my $r = ref($a) eq 'ARRAY' } } ],
  [ 'form(), current: ndarray',                    sub { for (1..$K) { my $f = form_now($p) } } ],
  [ 'form(), isa first: ndarray',                  sub { for (1..$K) { my $f = form_isa($p) } } ],
  [ 'form(), current: single value',               sub { for (1..$K) { my $f = form_now($a) } } ],
  [ 'form(), isa first: single value',             sub { for (1..$K) { my $f = form_isa($a) } } ],
  [ 'form(), current: array of values',            sub { for (1..$K) { my $f = form_now($aoa) } } ],
  [ 'form(), isa first: array of values',          sub { for (1..$K) { my $f = form_isa($aoa) } } ],
  [ 'form(), current: columns',                    sub { for (1..$K) { my $f = form_now($h) } } ],
  [ 'form(), isa first: columns',                  sub { for (1..$K) { my $f = form_isa($h) } } ],
);
die "forms disagree\n" if grep { form_now($_) != form_isa($_) } $p, $a, $aoa, $h;
printf "%-44s %8s\n", 'operation', 'ns';
printf "%-44s %8.1f\n", $_->[0], bench($_->[1]) for @rows;
