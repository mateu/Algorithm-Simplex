use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use PDL::Lite;
use Algorithm::Simplex::PDL;
use Data::Dumper;

=head1 Name

solve_PDL.pl - solve a PDL LP using custom solver.

=cut

# Get shell tableau object for access to EPSILON and MAXIMUM_PIVOTS
my $tableau_shell = Algorithm::Simplex->new( tableau => [ [] ] );

my $LP = {
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
  solve_LP( 'piddle', $LP->{'Baumol Advertising'}->{'initial_tableau'} );
print "Finished.\n";

=head1 Subroutines

=head2 solve_LP

Custom made solver for PDL.

=cut

sub solve_LP {
    my $model   = shift;
    my $tableau = shift;
    
    $tableau = pdl $tableau;

    my $tableau_object = Algorithm::Simplex::PDL->new(tableau => $tableau);

    my $counter = 1;
    until ( $tableau_object->is_optimal ) {
        my ( $pivot_row_number, $pivot_column_number ) =
          $tableau_object->determine_bland_pivot_row_and_column_numbers;
        $tableau_object->pivot( $pivot_row_number, $pivot_column_number );
        $tableau_object->exchange_pivot_variables( $pivot_row_number,
            $pivot_column_number );
        $counter++;
        die
"HALT: Exceeded the maximum number of pivots allowed: ". $tableau_shell->MAXIMUM_PIVOTS
          if ( $counter > $tableau_shell->MAXIMUM_PIVOTS );
    }
    
    return $tableau_object;
}
