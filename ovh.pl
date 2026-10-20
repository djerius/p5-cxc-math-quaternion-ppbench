use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# Per-call fixed costs (ns per call) at small N, to see where the difference
# between "kernel allocates the output" and "caller allocates it" comes from.

my $q = pdl(0.1,-0.3,0.5,0.8); $q /= sqrt(($q**2)->sum);
sub bench { my ($c)=@_; $c->() for 1..100; my ($n,$t0)=(0,time); do { $c->() for 1..100; $n+=100 } while time-$t0 < 1; (time-$t0)/$n }
my $D = double;

printf "%-44s %8s %8s\n", 'operation (ns per call)', 'N=1', 'N=1000';
my %v = map { $_ => random(3,$_) - 0.5 } 1, 1000;
my %r = map { $_ => zeroes(3,$_) } 1, 1000;
my @rows = (
  ['zeroes(3,N)                  [Core.pm wrapper]', sub { zeroes(3,$_[0]) }],
  ['PDL->zeroes(3,N)             [XS directly]',     sub { PDL->zeroes(3,$_[0]) }],
  ['PDL->new_from_specification($D,3,N)  [no data]', sub { PDL->new_from_specification($D,3,$_[0]) }],
  ['  ... ->make_physical',                          sub { PDL->new_from_specification($D,3,$_[0])->make_physical }],
  ['rotq($q,$v)                  [kernel allocs]',   sub { rotq($q,$v{$_[0]}) }],
  ['rotq($q,$v,$r)               [preallocated]',    sub { rotq($q,$v{$_[0]},$r{$_[0]}) }],
  ['rotq($q,$v,zeroes(3,N))',                        sub { rotq($q,$v{$_[0]},zeroes(3,$_[0])) }],
  ['rotq($q,$v,PDL->zeroes(3,N))',                   sub { rotq($q,$v{$_[0]},PDL->zeroes(3,$_[0])) }],
  ['rotq($q,$v,PDL->new_from_spec($D,3,N))',         sub { rotq($q,$v{$_[0]},PDL->new_from_specification($D,3,$_[0])) }],
);
for my $row (@rows) {
    my ($label, $code) = @$row;
    printf "%-44s %8.0f %8.0f\n", $label, map { my $N = $_; 1e9 * bench(sub { $code->($N) }) } 1, 1000;
}
