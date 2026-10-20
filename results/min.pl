# Per-cell minimum over several runs of one benchmark script.
#
#   perl results/min.pl results/cols-default-?.txt > results/cols-default-min.txt
#
# Lines are matched by position. Every plain decimal number (123, 4.5) is
# replaced by its minimum over the runs, keeping the first run's width and
# number of decimals; everything else (labels, N in 1e+03 form, headers)
# is copied from the first run.
use strict; use warnings;

my @runs = map { open my $h, '<', $_ or die "$_: $!\n"; [ <$h> ] } @ARGV;
die "usage: $0 run1.txt run2.txt ...\n" unless @runs;

my $num = qr/(?<![\w.+-])\d+(?:\.\d+)?(?![\w.+-])/;
for my $i (0 .. $#{ $runs[0] }) {
    my $line = $runs[0][$i];
    my @cols = map { [ ($_->[$i] // '') =~ /$num/g ] } @runs;
    my $k = 0;
    $line =~ s{$num}{
        my $tok = $&;
        my ($min) = sort { $a <=> $b } grep { defined } map { $_->[$k] } @cols;
        $k++;
        my $dp = $tok =~ /\.(\d+)/ ? length $1 : 0;
        sprintf '%*.*f', length $tok, $dp, $min;
    }ge;
    print $line;
}
