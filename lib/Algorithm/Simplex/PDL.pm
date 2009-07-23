package Algorithm::Simplex::PDL;
use Moose;
extends 'Algorithm::Simplex';
with 'Algorithm::Simplex::Role::Solve';
use Algorithm::Simplex::Types;
use PDL::Lite;
use namespace::autoclean;

=head1 Name

Algorithm::Simplex::PDL - PDL model of the Simplex Algorithm

=cut

# TODO: Probably need EPSILON for zero approximation check like in Float model.

has '+tableau' => (
    isa        => 'Piddle',
    coerce     => 1,
);

has '+display_tableau' => (
    isa        => 'PiddleDisplay',
    coerce     => 1,
);

=head1 Methods

=head2 _build_number_of_rows 

Set the number of rows.  This is actually for the A matrix in Ax <= y.
So the number is one less than the total number of rows in the tableau.
The same holds for number of columns.

=cut

sub _build_number_of_rows {
    my $self = shift;
    my ( $number_of_columns, $number_of_rows ) = ( $self->tableau->dims);
    return $number_of_rows - 1;
}

=head2 _build_number_of_columns 

set the number of columns given the tableau matrix

=cut

sub _build_number_of_columns {
    my $self = shift;
    my ( $number_of_columns, $number_of_rows ) = ( $self->tableau->dims );
    return $number_of_columns - 1;
}

=head2 pivot

Do the algebra of a Tucker/Bland pivot.  i.e. Traverse from one node to and 
adjacent node along the Simplex of feasible solutions.  This pivot method
is particular to this PDL model.

=cut

sub pivot {
    my $self                = shift;
    my $pivot_row_number    = shift;
    my $pivot_column_number = shift;

    my $pdl_A   = $self->tableau;
    my $neg_one = PDL->zeroes(1);
    $neg_one -= 1;

    my $scale_copy =
      $pdl_A->slice("($pivot_column_number),($pivot_row_number)")->copy;
    my $scale = $pdl_A->slice("($pivot_column_number),($pivot_row_number)");
    my $pivot_row = $pdl_A->slice(":,($pivot_row_number)");
    $pivot_row /= $scale_copy;
    $scale     /= $scale_copy;

    # peform pivot algebra in non-pivot rows
    for my $i ( 0 .. $self->number_of_rows ) {
        if ( $i != $pivot_row_number ) {
            my $a_ic_copy =
              $pdl_A->slice("($pivot_column_number),($i)")->copy;
            my $a_ic       = $pdl_A->slice("($pivot_column_number),($i)");
            my $change_row = $pdl_A->slice(":,($i)");
            my $diff_term  = $a_ic x $pivot_row;
            $change_row -= $diff_term;
            my $tmp = $neg_one x $a_ic_copy;
            $a_ic .= $tmp;    # $scale_copy;
            $a_ic /= $scale_copy;
        }
    }

    return $pdl_A;
}

# Count pivots made 
after 'pivot' => sub {
    my $self = shift;
    $self->number_of_pivots_made( $self->number_of_pivots_made + 1 );
    return;
};

=head2 is_optimal

Return 1 if the current solution is optimal, 0 otherwise.

=cut

sub is_optimal {
    my $self  = shift;
    my $T_pdl = $self->tableau;

    # Look at basement row to see if no positive entries exists.
    my $n_cols_A       = $self->number_of_columns - 1;
    my $number_of_rows = $self->number_of_rows;
    my $basement_row   = $T_pdl->slice("0:$n_cols_A,($number_of_rows)");
    my @basement_row   = $basement_row->list;
    my @positive_profit_column_numbers;
    foreach my $profit_coefficient (@basement_row) {
        if ( $profit_coefficient > 0 ) {
            return 0;
        }
    }

    return 1;
}

=head2 determine_simplex_pivot_columns

Look at the basement row to see where positive entries exists. 
Columns with positive entries in the basement row are pivot column candidates.

Should run optimality test, is_optimal, first to insure 
at least one positive entry exists in the basement row which then 
means we can increase the objective value for the maximization problem.

=cut 

sub determine_simplex_pivot_columns {
    my $self = shift;

    my @simplex_pivot_column_numbers;
    my $n_cols_A       = $self->number_of_columns - 1;
    my $number_of_rows = $self->number_of_rows;
    my $basement_row = $self->tableau->slice("0:$n_cols_A,($number_of_rows)");
    my @basement_row = $basement_row->list;
    my $column_number = 0;
    foreach my $profit_coefficient (@basement_row) {

        if ( $profit_coefficient > 0 ) {
            push @simplex_pivot_column_numbers, $column_number;
        }
        $column_number++;
    }

    return @simplex_pivot_column_numbers;
}

=head2 determine_positive_ratios

Starting with the pivot column find the entry that yields the lowest
positive b to entry ratio that has lowest bland number in the event of ties.

=cut

sub determine_positive_ratios {
    my $self                = shift;
    my $pivot_column_number = shift;

    my $n_rows_A          = $self->number_of_rows - 1;
    my $number_of_columns = $self->number_of_columns;
    my $pivot_column =
      $self->tableau->slice("($pivot_column_number),0:$n_rows_A");
    my @pivot_column = $pivot_column->list;
    my $constant_column =
      $self->tableau->slice("($number_of_columns),0:$n_rows_A");
    my @constant_column = $constant_column->list;
    my $row_number      = 0;
    my @positive_ratio_row_numbers;
    my @positive_ratios;

    foreach my $i ( 0 .. $n_rows_A ) {
        if ( $pivot_column[$i] > 0 ) {
            push @positive_ratios,
              ( $constant_column[$i] / $pivot_column[$i] );
            push @positive_ratio_row_numbers, $i;
        }
    }
    return ( \@positive_ratios, \@positive_ratio_row_numbers );
}

=head2 display_pdl

Given a Piddle return it as a string in a Matrix like format.

=cut

sub display_pdl {
    my $self   = shift;
    my $pdl    = $self->tableau;
    my $output = "$pdl";
    return $output;
}

=head2 current_solution

Return both the primal (max) and dual (min) solutions for the tableau.

=cut

sub current_solution {
    my $self = shift;

    # Report the Current Solution as primal dependents and dual dependents.
    my @y = @{ $self->y_variables };
    my @u = @{ $self->u_variables };

    # Dependent Primal Variables
    my $n_rows_A          = $self->number_of_rows - 1;
    my $number_of_columns = $self->number_of_columns;
    my $constant_column =
      $self->tableau->slice("($number_of_columns),0:$n_rows_A");
    my @constant_column = $constant_column->list;
    my %primal_solution;
    for my $i ( 0 .. $#y ) {
        $primal_solution{ $y[$i]->{generic} } = $constant_column[$i];
    }

    # Dependent Dual Variables
    my $n_cols_A       = $self->number_of_columns - 1;
    my $number_of_rows = $self->number_of_rows;
    my $basement_row = $self->tableau->slice("0:$n_cols_A,($number_of_rows)");
    my @basement_row = $basement_row->list;
    my %dual_solution;
    for my $j ( 0 .. $#u ) {
        $dual_solution{ $u[$j]->{generic} } = $basement_row[$j] * (-1);
    }

    return ( \%primal_solution, \%dual_solution );
}

__PACKAGE__->meta->make_immutable;

1;
