use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/blib/lib", "$FindBin::Bin/blib/arch", $FindBin::Bin;
use PDL; use QBench; use QDispatch; use Scalar::Util qw(refaddr);

# Correctness of the two dispatch prototypes against the kernels they call.
my @q = (0.1, -0.3, 0.5, 0.8); { my $n = sqrt(eval { my $s = 0; $s += $_**2 for @q; $s }); $_ /= $n for @q }
my ($fail, $n) = (0, 0);
sub ok { my ($c, $name) = @_; $n++; unless ($c) { $fail++; print "FAIL: $name\n" } }
sub close_enough { my ($a, $b) = @_; (abs(pdl($a) - pdl($b)))->max < 1e-12 }

my $vp  = random(3, 7) - 0.5;
my $aoa = unpdl($vp);
my %col = ( x => [ $vp->slice('(0)')->list ], y => [ $vp->slice('(1)')->list ], z => [ $vp->slice('(2)')->list ] );
my $ref = unpdl(rotq(pdl(@q), $vp));
my $col_ref = { x => [ map { $_->[0] } @$ref ], y => [ map { $_->[1] } @$ref ], z => [ map { $_->[2] } @$ref ] };
my $col_ok = sub { my ($h) = @_; close_enough([ @{$h}{qw(x y z)} ], [ @{$col_ref}{qw(x y z)} ]) && @{ $h->{x} } == 7 };

for my $d ([ perl => \&QDispatch::rotate_vec ], [ xs => \&QBench::xsd_rotate_vec ]) {
    my ($name, $f) = @$d;
    # single vector
    ok(close_enough($f->(\@q, $aoa->[0]), $ref->[0]), "$name single new");
    my @o = (9, 9, 9, 9, 9); my $r = $f->(\@q, $aoa->[0], \@o);
    ok(refaddr($r) == refaddr(\@o) && @o == 3 && close_enough(\@o, $ref->[0]), "$name single into, truncated");
    my @s = @{ $aoa->[0] }; $f->(\@q, \@s, \@s); ok(close_enough(\@s, $ref->[0]), "$name single in place");
    # records
    ok(close_enough($f->(\@q, $aoa), $ref), "$name records new");
    my @ro = map { [ 1, 2, 3 ] } 1..9; $f->(\@q, $aoa, \@ro);
    ok(@ro == 7 && close_enough(\@ro, $ref), "$name records into, outer truncated");
    if ($name eq 'xs') {   # the Perl prototype calls the old kernel, which does not shorten inner records
        my @r4 = map { [ 1, 2, 3, 4 ] } 1..7; $f->(\@q, $aoa, \@r4);
        ok(!grep({ @$_ != 3 } @r4) && close_enough(\@r4, $ref), "$name records into, inner records truncated");
    }
    my $w = [ map { [@$_] } @$aoa ]; $f->(\@q, $w, $w); ok(close_enough($w, $ref), "$name records in place");
    # columns
    ok($col_ok->($f->(\@q, \%col)), "$name columns new");
    my %co = ( x => [ (0) x 9 ], extra => 1 ); $f->(\@q, \%col, \%co);
    ok($col_ok->(\%co) && $co{extra} == 1, "$name columns into, truncated, other keys kept");
    my %cw = map { $_ => [ @{ $col{$_} } ] } qw(x y z); $f->(\@q, \%cw, \%cw); ok($col_ok->(\%cw), "$name columns in place");
    # ndarrays (the XS prototype hands these to the Perl dispatcher)
    ok(close_enough(unpdl($f->(\@q, $vp)), $ref), "$name ndarray new");
    my $po = zeroes(3, 7); $f->(\@q, $vp, $po); ok(close_enough(unpdl($po), $ref), "$name ndarray into");
    my $pw = $vp->copy; $f->(\@q, $pw, $pw); ok(close_enough(unpdl($pw), $ref), "$name ndarray in place");
    ok(!eval { $f->(\@q, $vp, zeroes(float, 3, 7)); 1 } && $@ =~ /double/, "$name ndarray float output croaks");
    # errors
    ok(!eval { $f->(\@q, $aoa, {}); 1 }, "$name records into a hash croaks");
    ok(!eval { $f->(\@q, [ [ 1, 2 ] ]); 1 }, "$name short record croaks");
    ok(!eval { $f->(\@q, \"x"); 1 }, "$name non-container croaks");
}
printf "%d/%d checks passed\n", $n - $fail, $n;
