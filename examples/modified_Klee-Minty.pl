use Algorithm::Simplex::Rational;
use Data::Dumper;

my $klee_minty = [
    [ 1,  0, 0, 0, 1 ],
    [ 4,  1, 0, 0, 8 ],
    [ 8,  4, 1, 0, 64 ],
    [ 16, 8, 4, 1, 512 ],
    [ 8,  4, 2, 1, 0 ]
];

#my $lumber_mill = [ [ 1, 3, 2, 10 ], [ 2, 1, 1, 8 ], [ 3, 2, 4, 0 ] ];
#my $ad = [ [ 8, 3, 4, 40 ], [ 40, 10, 10, 200 ], [ 160, 60, 80, 0 ], ];
#my $bland = [
#    [ 1 / 4, -8,  -1,     9,  0 ],
#    [ 1 / 2, -12, -1 / 2, 3,  0 ],
#    [ 0,     0,   1,      0,  1 ],
#    [ 3 / 4, -20, 1 / 2,  -6, 0 ],
#];

my $problem = Algorithm::Simplex::Rational->new( tableau => $klee_minty );
if ( $problem->solve ) {
    print Dumper $problem->current_solution;
    print Dumper $problem->display_tableau;
    print 'Objective function value: ', $problem->objective_function_value, "\n";
    print 'Number of pivots made: ', $problem->number_of_pivots_made, "\n";
}
else {
    print "Bitter.\n";
}
