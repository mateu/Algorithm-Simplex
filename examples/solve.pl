use FindBin;
use lib "$FindBin::Bin/../lib";
use Algorithm::Simplex::Rational;
use Algorithm::Simplex::PDL;
use Algorithm::Simplex::Float;
use PDL::Lite;
use Data::Dumper;

=head1 Name

solve.pl - a script that demonstrates the solve() method on each model.

=cut

my $LP = {
    'McRae: Lumber Mill' => {
        initial_tableau =>
          [ [ 1, 3, 2, 10 ], [ 2, 1, 1, 8 ], [ 3, 2, 4, 0 ] ],
        optimal_tableau => [
            [ -1 / 3, 2 / 3,  5 / 3,   4 ],
            [ 2 / 3,  -1 / 3, -1 / 3,  2 ],
            [ -2 / 3, -5 / 3, -11 / 3, -22 ]
        ],
    },
};

my $initial_tableau = $LP->{'McRae: Lumber Mill'}->{'initial_tableau'};

# Copy by value
my $initial_tableau_2 = matrix_copy($initial_tableau);
my $initial_tableau_3 = matrix_copy($initial_tableau);
my $tableau_object;

$tableau_object = Algorithm::Simplex::PDL->new( tableau => $initial_tableau );
if ( my $final_tableau_object = $tableau_object->solve ) {
    print Dumper $final_tableau_object->current_solution;
}
else {
    print "Exceeding maximum number of allowed loops: "
      . $tableau_object->MAXIMUM_PIVOTS, "\n";
}

$tableau_object =
  Algorithm::Simplex::Rational->new( tableau => $initial_tableau_2 );
if ( my $final_tableau_object = $tableau_object->solve ) {
    print Dumper $final_tableau_object->current_solution;
}
else {
    print "Exceeding maximum number of allowed loops: "
      . $tableau_object->MAXIMUM_PIVOTS, "\n";
}

$tableau_object =
  Algorithm::Simplex::Float->new( tableau => $initial_tableau_3 );
if ( my $final_tableau_object = $tableau_object->solve ) {
    print Dumper $final_tableau_object->current_solution;
}
else {
    print "Exceeding maximum number of allowed loops: "
      . $tableau_object->MAXIMUM_PIVOTS, "\n";
}

sub matrix_copy {

    # code to copy matrix
    my $matrix = shift;
    my $matrix_copy;
    my $rows    = scalar @{$matrix};
    my $columns = scalar @{ $matrix->[0] };

    for my $i ( 0 .. $rows - 1 ) {
        for my $j ( 0 .. $columns - 1 ) {
            $matrix_copy->[$i]->[$j] = $matrix->[$i]->[$j];
        }
    }

    return $matrix_copy;
}

