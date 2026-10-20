# per-cell minimum over several runs of cols.pl output; non-numeric lines from the first file
my @f = map { [ do { open my $h, '<', $_ or die; <$h> } ] } @ARGV;
for my $i (0 .. $#{$f[0]}) {
    my $l = $f[0][$i];
    if ($l =~ /^\s+\de[+-]\d\d\s/) {
        my @rows = map { [ split ' ', $_->[$i] ] } @f;
        my @m = map { my $c = $_; my ($min) = sort { $a <=> $b } map { $_->[$c] } @rows; $min } 1 .. $#{$rows[0]};
        my @w = map { length } ($l =~ /(\s+\S+)/g);
        printf "%*s" x @w . "\n", map { ($w[$_], !$_ ? $rows[0][0] : $rows[0][$_] =~ /^[\d.]+$/ ? sprintf('%.1f',$m[$_-1]) : $rows[0][$_]) } 0..$#w;
    } else { print $l }
}
