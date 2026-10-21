use strict; use warnings;

# Evaluates: "multiplying N quaternions by 1 vector works if you vertically
# stack the quaternion (rotation) matrices, then do a splitdim on the output
# extremely tall vector".  Checks it against rotq, and times it against
# broadcasting the matrices and against rotq broadcasting the quaternions.
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch";
use PDL; use QBench; use Time::HiRes qw(time);

# rotation matrices (3,3,N) from quaternions (4,N), PDL layout: $m->at(col, row, i)
sub matrices {
    my ($q) = @_;
    my ($x, $y, $z, $w) = map { $q->slice("($_)") } 0..3;
    my $m = zeroes(3, 3, $q->dim(1));
    $m->slice('(0),(0)') .= 1 - 2*($y*$y + $z*$z); $m->slice('(1),(0)') .= 2*($x*$y - $z*$w); $m->slice('(2),(0)') .= 2*($x*$z + $y*$w);
    $m->slice('(0),(1)') .= 2*($x*$y + $z*$w);     $m->slice('(1),(1)') .= 1 - 2*($x*$x + $z*$z); $m->slice('(2),(1)') .= 2*($y*$z - $x*$w);
    $m->slice('(0),(2)') .= 2*($x*$z - $y*$w);     $m->slice('(1),(2)') .= 2*($y*$z + $x*$w);     $m->slice('(2),(2)') .= 1 - 2*($x*$x + $y*$y);
    $m;
}
sub bench { my ($c) = @_; $c->(); my ($n, $t0) = (0, time); do { $c->(); $n++ } while time - $t0 < 0.5; 1e9 * (time - $t0) / $n }

my $N = 5;
my $q = random(4, $N) - 0.5; $q /= sqrt(($q**2)->sumover)->dummy(0);
my $v = pdl(0.3, -0.2, 0.9);
my $ref = rotq($q, $v);                                  # (3,N): q broadcast against one vector

my $m     = matrices($q);                                # (3,3,N)
my $stack = $m->clump(1, 2);                             # (3,3N): vertical stack, a view
my $tall  = $stack x $v->dummy(0);                       # (3,3N) x (1,3) column -> (1,3N)
my $split = $tall->slice('(0)')->splitdim(0, 3);         # (3N) -> (3,N)
printf "stack dims %s, product dims %s, split dims %s\n", join(',', $stack->dims), join(',', $tall->dims), join(',', $split->dims);
printf "stack+splitdim matches rotq:   max |diff| = %.1e\n", (abs($split - $ref))->max;
my $bc = ($m x $v->dummy(0))->slice('(0)');              # (3,3,N) x (1,3) broadcasts -> (1,3,N)
printf "broadcast matrices matches:    max |diff| = %.1e   dims %s\n", (abs($bc - $ref))->max, join(',', $bc->dims);
printf "clump is a view (no copy):     %s\n", $stack->is_inplace ? '?' : ($stack->get_dataref == $m->get_dataref ? 'shares data' : 'virtual');

print "\nns per quaternion (one vector, N quaternions):\n";
printf "%8s %10s %14s %14s %16s\n", 'N', 'rotq', 'matrices only', 'stack+x+split', 'broadcast x';
for my $n (1e2, 1e4, 1e6) {
    my $qq = random(4, $n) - 0.5; $qq /= sqrt(($qq**2)->sumover)->dummy(0);
    my $mm = matrices($qq);
    my @t = map { bench($_) / $n }
        sub { rotq($qq, $v) },
        sub { matrices($qq) },
        sub { ($mm->clump(1, 2) x $v->dummy(0))->slice('(0)')->splitdim(0, 3)->sever },
        sub { ($mm x $v->dummy(0))->slice('(0)')->sever };
    printf "%8.0e %10.1f %14.1f %14.1f %16.1f\n", $n, @t;
}
