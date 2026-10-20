use strict; use warnings;
use PDL; use Time::HiRes qw(time);
# Split the ~1 us cost of creating a Perl-visible ndarray into construction vs. destruction.
my $D = double; my $K = 200_000;
sub t { my ($c) = @_; my $t0 = time; $c->(); 1e9 * (time - $t0) / $K }
for my $rep (1..3) {
    my @keep;
    my $create  = t(sub { push @keep, PDL->new_from_specification($D,3,1) for 1..$K });
    my $destroy = t(sub { @keep = () });
    my $null    = t(sub { my $x; $x = PDL->null for 1..$K });
    my $nfsnt   = t(sub { my $x; $x = PDL->new_from_specification(3,1) for 1..$K });
    my $sub     = t(sub { my $x; $x = main->can('t') for 1..$K });
    printf "create %4.0f  destroy %4.0f  | null create+destroy %4.0f  nfs(no type) c+d %4.0f  (loop/method baseline %3.0f) ns\n",
        $create, $destroy, $null, $nfsnt, $sub;
}
