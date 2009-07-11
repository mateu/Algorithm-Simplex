package Algorithm::Simplex::Role::Solve;
use Moose::Role;

=head1 Name

Algorithm::Simplex::Role::Solve - solve() method implemented as Moose role.

=cut

=head1 Synposis

    use Algorithm::Simplex::Rational;
    use Data::Dumper;
    my $matrix = [
        [ 5,  2,  30],
        [ 3,  4,  20],
        [10,  8,   0],
    ];
    my $tableau_object = Algorithm::Simplex::Rational->new( tableau => $matrix );
    $tableau_object->solve;
    print Dumper $tableau_object->display_tableau;
    my ($primal_solution, $dual_solution) = $tableau_object->current_solution;
    print Dumper $primal_solution;
    print Dumper $dual_solution;
     
=cut    

requires 'tableau', 
         'determine_bland_pivot_row_and_column_numbers',
         'pivot',
         'exchange_pivot_variables';

sub solve {
    my $tableau_object = shift;

    my $counter = 1;
    until ( $tableau_object->is_optimal ) {
        my ( $pivot_row_number, $pivot_column_number ) =
          $tableau_object->determine_bland_pivot_row_and_column_numbers;
        $tableau_object->pivot( $pivot_row_number, $pivot_column_number );
        $tableau_object->exchange_pivot_variables( $pivot_row_number,
            $pivot_column_number );
        $counter++;

        # Too many pivots?
        if ( $counter > $tableau_object->MAXIMUM_PIVOTS ) {
            warn "HALT: Exceeded the maximum number of pivots allowed: "
              . $tableau_object->MAXIMUM_PIVOTS . "\n";
            return 0;
        }
    }

    return 1;
}

1;
