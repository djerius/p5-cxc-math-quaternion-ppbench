use PDL::Core::Dev; use ExtUtils::MakeMaker;
my $pkg = ['qbench.pd', 'QBench', 'QBench', undef, 1];
WriteMakefile(pdlpp_stdargs($pkg), OPTIMIZE => '-O2');
sub MY::postamble { pdlpp_postamble($pkg) }
