package Algorithm::Simplex::PDL;
use Moose;
use Algorithm::Simplex::Types;
use namespace::autoclean;
extends 'Algorithm::Simplex';
with 'Algorithm::Simplex::Role::Solve';
use PDL::Lite;



=head1 Name

Algorithm::Simplex::PDL - PDL model of the Simplex Algorithm

=head1 Methods

=head2 pivot

Do the algebra of a Tucker/Bland pivot.  i.e. Traverse from one node to and 
adjacent node along the Simplex of feasible solutions.

=cut

# TODO: Probably need EPSILON for zero approximation check like in Float model.

has tableau => (
    is       => 'rw',
    isa      => 'Piddle',
    required => 1,
    coerce   => 1,
);

=head1 Methods

=head2 _build_number_of_rows 

Set the number of rows.  This is actually for the A matrix in Ax <= y.
So the number is one less than the total number of rows in the tableau.
The same holds for number of columns.

=cut

sub _build_number_of_rows {
    my $self = shift;
    my ( $number_of_columns, $number_of_rows ) = PDL::dims( $self->tableau );
    return $number_of_rows - 1;
}

=head2 _build_number_of_columns 

set the number of columns given the tableau matrix

=cut

sub _build_number_of_columns {
    my $self = shift;
    my ( $number_of_columns, $number_of_rows ) = PDL::dims( $self->tableau );
    return $number_of_columns - 1;
}

no Moose;

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

sub is_optimal {
    my $self  = shift;
    my $T_pdl = $self->tableau;

    # Look at basement row to see if no positive entries exists.
    my $n_cols_A       = $self->number_of_columns - 1;
    my $number_of_rows = $self->number_of_rows;
    my $basement_row   = $T_pdl->slice("0:$n_cols_A,($number_of_rows)");
    my @basement_row   = $basement_row->list;
    my @positive_profit_column_numbers;
    my $optimal_flag = 1;
    foreach my $profit_coefficient (@basement_row) {
        if ( $profit_coefficient > 0 ) {
            $optimal_flag = 0;
            last;
        }
    }

    return $optimal_flag;
}

sub determine_simplex_pivot_columns {
    my $self = shift;

    my $T_pdl = $self->tableau;

    # Look at basement row to see where positive entries exists.
    # Run optimality test first to insure at least one positve profit exists.
    my @simplex_pivot_column_numbers;
    my $n_cols_A       = $self->number_of_columns - 1;
    my $number_of_rows = $self->number_of_rows;
    my $basement_row   = $T_pdl->slice("0:$n_cols_A,($number_of_rows)");
    my @basement_row   = $basement_row->list;
    my $column_number  = 0;
    foreach my $profit_coefficient (@basement_row) {

        if ( $profit_coefficient > 0 ) {
            push @simplex_pivot_column_numbers, $column_number;
        }
        $column_number++;
    }

    return @simplex_pivot_column_numbers;
}

sub determine_positive_ratios {

# Starting with the pivot column find the entry that yields the lowest
# positive b to entry ratio that has lowest bland number in the event of ties.
    my $self                = shift;
    my $pivot_column_number = shift;

    my $pdl               = $self->tableau;
    my $n_rows_A          = $self->number_of_rows - 1;
    my $number_of_columns = $self->number_of_columns;
    my $pivot_column      = $pdl->slice("($pivot_column_number),0:$n_rows_A");
    my @pivot_column      = $pivot_column->list;
    my $constant_column   = $pdl->slice("($number_of_columns),0:$n_rows_A");
    my @constant_column   = $constant_column->list;
    my $row_number        = 0;
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

sub get_pdl {
    my $self   = shift;
    my $pdl    = $self->tableau;
    my $output = "$pdl";
    return $output;
}

__PACKAGE__->meta->make_immutable;

1;
