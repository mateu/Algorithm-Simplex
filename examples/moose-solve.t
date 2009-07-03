use Algorithm::Simplex::Float;
use Data::Dumper;

$MAXIMUM_PIVOTS = 200;

$LP = {
    'Baumol Advertising' => {
        initial_tableau =>
          [ [ 8, 3, 4, 40 ], [ 40, 10, 10, 200 ], [ 160, 60, 80, 0 ], ],
        optimal_tableau => [
            [ 1 / 8, 3 / 8, 1 / 2, 5 ],
            [ -5,    -5,    -10,   0 ],
            [ -20,   0,     0,     -800 ],
        ],
    },
};

my $final_tableau_object =
  solve_LP( 'float', $LP->{'Baumol Advertising'}->{'initial_tableau'} );

print "Finished.\n";
#print Dumper $final_tableau_object->tableau;  


sub solve_LP {
    my $model   = shift;
    my $tableau = shift;

    my $tableau_object =
      $model eq 'float'
      ? Algorithm::Simplex::Float->new(tableau => $tableau)
      : die "The model type: $model could not be found.";

    my $counter = 1;
    until ( $tableau_object->is_optimal ) {
        my ( $pivot_row_number, $pivot_column_number ) =
          $tableau_object->determine_bland_pivot_row_and_column_numbers;
        $tableau_object->pivot( $pivot_row_number, $pivot_column_number );
        $tableau_object->exchange_pivot_variables( $pivot_row_number,
            $pivot_column_number );
        $counter++;
        die
"HALT: Exceeded the maximum number of pivots allowed: ". $tableau_object->MAXIMUM_PIVOTS
          if ( $counter > $tableau_object->MAXIMUM_PIVOTS );
    }
    
    return $tableau_object;
}
