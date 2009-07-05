package Algorithm::Simplex::Role::Solve;
use Moose::Role;

=head1 Name

Algorithm::Simplex::Role::Solve - solve() method implemented as Moose role.

=cut

=head1 Synposis

    use Algorithm::Simplex::Float;
    use Data::Dumper;
    my $matrix =  [ [ 1, 3, 2, 10 ], [ 2, 1, 1, 8 ], [ 3, 2, 4, 0 ] ];
    my $tableau = Algorithm::Simplex::Float->new( tableau => $matrix );
    my $final_tableau = $tableau->solve;
    print Dumper $final_tableau->tableau . "\n";
    
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

    return $tableau_object;
}

1;
