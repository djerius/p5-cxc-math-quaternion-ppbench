# One case for instruction counting under callgrind:
#   perl ir-case.pl CASE VARIANT K
# Runs the case K times after setting up its data, so that the counts for
# K and for 0 calls differ only by the K calls.
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench;
my ($case, $variant, $K) = @ARGV;
my @q = (0.1, -0.3, 0.5, 0.8); { my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my $N = 1000;
my $vp = sequence(3, $N) / (3 * $N) - 0.5; my $aoa = unpdl($vp);
my %col = map { my $i = $_; ( (qw(x y z))[$i] => [ $vp->slice("($i)")->list ] ) } 0..2;
my @v = (0.3, -0.2, 0.9); my @o = (0, 0, 0); my @ro = map { [0, 0, 0] } 1..$N; my %co = map { $_ => [ (0) x $N ] } qw(x y z);
my $w = [ map { [@$_] } @$aoa ]; my %cw = map { $_ => [ @{ $col{$_} } ] } qw(x y z);
my %run = (
  'single into'      => { kernel => sub { xs_rot1_into(\@q, \@v, \@o) }, xsd => sub { xsd_rotate_vec(\@q, \@v, \@o) }, qd => sub { qd_rotate_vec(\@q, \@v, \@o) } },
  'single new'       => { kernel => sub { my $r = xs_rot1(\@q, \@v) }, xsd => sub { my $r = xsd_rotate_vec(\@q, \@v) }, qd => sub { my $r = qd_rotate_vec(\@q, \@v) } },
  '1000 values into' => { kernel => sub { xs_rot_aoa_into(\@q, $aoa, \@ro) }, xsd => sub { xsd_rotate_vec(\@q, $aoa, \@ro) }, qd => sub { qd_rotate_vec(\@q, $aoa, \@ro) } },
  '1000 values in place' => { kernel => sub { xs_rot_aoa_into(\@q, $w, $w) }, xsd => sub { xsd_rotate_vec(\@q, $w, $w) }, qd => sub { qd_rotate_vec(\@q, $w, $w) } },
  '1000 values new'  => { kernel => sub { my $r = xs_rot_aoa(\@q, $aoa) }, xsd => sub { my $r = xsd_rotate_vec(\@q, $aoa) }, qd => sub { my $r = qd_rotate_vec(\@q, $aoa) } },
  '1000 columns into' => { kernel => sub { xs_rot_cols_into(\@q, @col{qw(x y z)}, @co{qw(x y z)}) }, xsd => sub { xsd_rotate_vec(\@q, \%col, \%co) }, qd => sub { qd_rotate_vec(\@q, \%col, \%co) } },
  '1000 columns new' => { kernel => sub { my @r = xs_rot_cols(\@q, @col{qw(x y z)}) }, xsd => sub { my $r = xsd_rotate_vec(\@q, \%col) }, qd => sub { my $r = qd_rotate_vec(\@q, \%col) } },
);
my $f = $run{$case}{$variant} or die "no case '$case' / '$variant'\n";
$f->() for 1 .. $K;
